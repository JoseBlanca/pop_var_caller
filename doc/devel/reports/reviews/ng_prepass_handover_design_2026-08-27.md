# Review — the pre-pass → calling handover (`ng-prepass-handover`, 629e84ff..346abaaf)

Reading review, no code run. Everything below was checked against the tree at
`/Users/jose/devel/pop_var_caller-review-design` (detached at 346abaaf).

**Verdict.** The code does what the plan asked and the two admitted divergences are both
justified — I checked B1's premise myself and the report is right that it was false. The defects
are all in the prose: five separate factual claims about how the code behaves are wrong, and four
of them are the *stated reason* for a design or a test. One of them is contradicted by a sentence
three lines further down the same doc comment.

---

## Blocker 1 — "a missing minted-error total is not refused; the run finishes normally" is the opposite of what `assemble` does

**Where** (five places, same claim):

- report, step *The reads' own claim about their errors now leaves the fit*: "A read group with a
  fitted rate and no minted total is not refused downstream. It takes the defaulted calibration —
  scale one, every read of the library charged the error floor — and the run finishes normally with
  quietly overconfident reads."
- `PROJECT_STATUS.md`: "A library with a fitted error rate and no minted-error total is not
  refused: it takes the defaulted calibration, scale one, every read charged the error floor, and
  the run finishes normally with quietly overconfident reads."
- commit `1eafa835`: "A read group with a fitted rate and no minted total takes the defaulted
  calibration: scale one, every read charged the error floor, and the run completing normally."
- `src/ng/parameter_estimation/generic/mod.rs:457-461` (the new `minted_errors` doc comment).
- `src/ng/parameter_estimation/generic/estimate.rs:640-645` (the new test's doc comment).

**Why it is wrong.** `RunParameters::assemble` refuses exactly this case, and refused it *before*
this branch — `checked_read_group_count_of`, `src/ng/calling/run_parameters.rs:631-642`, present
unchanged at the branch point (`git show 629e84ff:src/ng/calling/run_parameters.rs`, line 458):

```rust
assert_eq!(
    error_rate_by_read_group.contains_key(group),
    minted_by_read_group.contains_key(group),
    "read group {} has {} — the fitted rate and the accumulator's total come from one \
     pass over one set of reads, so one without the other means they saw different data",
```

An empty `minted_errors` map with a non-empty `error_rate` map makes that assertion fail on the
first read group. The run does not finish; it aborts at assembly with a message naming the read
group. The same doc comment says so correctly four lines later
(`generic/mod.rs:466`: "The calling step's own assembly refuses a fitted rate whose read group has
no total here, and the reverse"), so the block contradicts itself.

**Second error inside the same sentence.** *Scale one* and *the error floor* are two different,
mutually exclusive outcomes. `ReadGroupCalibration::defaulted()`
(`src/ng/calling/likelihood/mod.rs:361-366`) is `scale: 1.0`, documented as "the qualities are used
as reported" — every read charged exactly what the instrument claimed. Charging every read the
floor is what a scale of **zero** does, which `from_fitted_rate`
(`likelihood/mod.rs:307-311`) exists to refuse.

**Replacement** (report and status; adapt for the commit and the two doc comments):

> Without this field nothing that assembled a run's parameters could supply the denominator, and
> `RunParameters::assemble` refuses that rather than absorbing it: `checked_read_group_count_of`
> asserts that a read group has a fitted rate and a minted total or neither, so a run whose
> per-sample results carried no totals aborts at assembly naming the first read group. The test is
> written against the quieter failure the assertion cannot see — a map filled from the wrong
> sample, whose keys are a different library's — so *the map is non-empty* is not what it asserts.

**Note for the doc comment at `generic/mod.rs:457`**: "The accumulators are dropped when this value
is returned" is also not true of both entry points. `GenericAccumulators::estimate` takes `&self`
(`estimate.rs:158`) and the new test in the same commit reads `driven.minted_errors()` *after* the
fit. Only the streaming entry point drops its local accumulator. Write: "The streaming entry point
drops its accumulator as soon as the fit returns, so this is the last moment the totals can be
reached."

---

## Blocker 2 — "a permuted per-sample list gives every sample a coefficient and a contamination fraction, just not its own"

**Where**: `src/ng/calling/run_parameters.rs:152-156` (the `from_prepass` contract), the report's
*Two joins that no type enforces* bullet, commit `d8970944` ("Nothing about that shows up as a
failure downstream"), commit `9ac41520`. It originates in the plan's D1 and is repeated as fact.

**Why it is wrong — two independent reasons.**

1. **Contamination is not read from the permuted list at all.** `from_prepass` looks it up by the
   *run's own* sample name, taken from `read_groups.read_groups_per_sample()`:
   `let of_the_fit = joint.contamination.get(sample)` (`run_parameters.rs:239`). `sample` comes
   from `of_sample.sample`, not from `generic_by_sample[index]`. Permuting the two per-sample lists
   cannot move a single contamination fraction.
2. **The coefficient misassignment is not silent — the seam panics on it.** In the same loop
   iteration, every key of `generic.error_rate`, `generic.minted_errors` and
   `repeat_tract.substitution_rate_by_stratum()` is passed through `its_own_read_group`. Read groups
   are disjoint across samples (a header's `@RG` names one `SM`), so a permuted list hands sample *i*
   a map keyed by somebody else's libraries and the assertion fires, naming the sample and the read
   group. The implementer's own defect table confirms this: "the samples read in reverse order — the
   tracing test, and the refusal test". The only surviving hole is a sample carrying **no**
   read-group-keyed values at all, and the contract does not mention that.

**Replacement** for the contract bullet:

> - **Sample order.** `generic_by_sample` and `repeat_tract_by_sample` are joined to the run's
>   sample table by position, while the cohort fit keys its own per-sample results by name — so the
>   contamination fractions follow the run's sample names and the coefficients follow the lists'
>   order, and only the second can be permuted. A permuted list would give a sample another
>   sample's coefficient; `its_own_read_group` catches it on that sample's very first error rate,
>   because read groups are disjoint across samples, and the run stops naming both. What it cannot
>   catch is a sample carrying no read-group-keyed value at all, which no pre-pass output produces
>   today.

---

## Major 3 — "it computes nothing: every number it hands on was measured by one of the three"

**Where**: report, *What it is, in one paragraph*; `PROJECT_STATUS.md` ("One function now does it
and computes nothing of its own"); `run_parameters.rs:141-143` ("**no number of its own**: every
value it hands on was computed by one of the three"); commit `9ac41520`.

**Why it is wrong.** `from_prepass` calls `Self::seed_from_moments`, which calls
`seed_from_population_moments` (`src/ng/calling/genotype_prior/seed_generic.rs:250`). That function
solves a concentration pair from the two moments (`total_for_diversity`) and, where a moment is
missing, substitutes constants the pre-pass never measured: `NEUTRAL_ALPHA_REF`,
`ExpectedHeterozygosity::SPECIES_FALLBACK`, `MIN_ALT_CONCENTRATION`. The module's own doc says so
plainly at `run_parameters.rs:118-120`: `seed_from_moments` is "the one number here that is
*derived* rather than gathered".

**Replacement**:

> It gathers rather than fits: every per-library and per-sample number it hands on is one of the
> three pre-pass outputs, unchanged. The one derived value is the genotype prior's seed, which
> `seed_from_moments` solves in closed form from the cohort fit's two moments — and falls back to
> its own constants where a moment is missing.

---

## Major 4 — "the constructor had ten callers and all ten were in its own test module"

**Where**: report, *What it is, in one paragraph*; `PROJECT_STATUS.md` ("its constructor had ten
callers and all ten were its own tests"). Inherited from the plan's opening paragraph.

**Why it is wrong.** At the branch point `RunParameters::assemble` has **29** call sites, in 29
distinct enclosing functions, plus a local `assemble(...)` helper called three more times:

```
git show 629e84ff:src/ng/calling/run_parameters.rs | grep -c "RunParameters::assemble("   # 29
```

The *interesting* half of the claim — that all of them are inside `mod tests` — is true.

**Replacement**: "the constructor had twenty-nine callers and every one of them was inside its own
test module, each handing it values written by hand."

---

## Major 5 — "the permutation left all 201 other tests in that file green"

**Where**: report, step *A sample's tract substitution rates*; commit `ebcd86a8` ("all 201 other
tests in the file green").

**Why it is wrong.** `src/ng/parameter_estimation/ssr/mod.rs` holds **91** `#[test]` functions, not
202. The number 203 is the whole `ssr` module directory (mod.rs 91, locus_offsets.rs 40,
stratum_table.rs 37, slippage.rs 29, offset_bucket.rs 6), so the count belongs to a `cargo test
ssr::` filter and not to a file — and it is 202 others, not 201.

**Replacement**: "the permutation left the other 202 tests of the `ssr` module green — the file's own
90 among them."

---

## Minor 6 — the thousand-sample table was measured through `assemble`, on a cohort with no contamination, and the report does not say so

**Where**: report, *One sample and a thousand*: "**A thousand** is a measurement, in
`examples/ng_prepass_handover_footprint.rs`, taken with dhat's allocator as live bytes."

Two things a reader would assume and should not:

- The example calls `RunParameters::assemble` and rebuilds the union by hand
  (`examples/ng_prepass_handover_footprint.rs:~265-300`), not `from_prepass`. The example's own
  comment says so; the report does not.
- It passes `&BTreeMap::new()` for contamination, so **the run-wide contamination map and the dense
  `contamination_by_read_group` vector are not in any column**. That is why "what assembling adds"
  is 24 bytes at a thousand samples: with every sample's fraction identified, assembling instead
  builds a `Vec<ContaminationView>` of one entry a library (32 bytes each, ≈32 kB at a thousand),
  and the union grows by a 1,000-entry `BTreeMap` of `ContaminationEstimate`. Small against 26.6 MB,
  but the "24 bytes" column reads as a property of the seam and is a property of an uncontaminated
  fixture.

**Replacement** (report): "**A thousand** is a measurement in
`examples/ng_prepass_handover_footprint.rs`, live bytes from dhat's allocator. It builds the union
the way `from_prepass` builds it and calls `assemble` directly, because `ReadGroups` can only be
minted from an alignment header or from the crate's own test-only constructor. Contamination is
left empty, so the last column — the 24 bytes assembling adds — is a run with no fractions
identified; a cohort where every library has one adds a dense view per library besides, 32 kB at a
thousand."

Also in the example's own header comment: "the run's read-group table … is read from the alignment
files' headers and cannot be built in memory" is not the reason. It *can* —
`ReadGroups::of_libraries` (`src/ng/read/input/read_groups.rs:423`) does exactly that — but it is
`#[cfg(test)] pub(crate)`, so an example cannot reach it. Say that instead.

---

## Minor 7 — "43 parts in 44"

**Where**: report and commit `2f35e8a4`: "never holds the other 43 parts in 44."

26,595,880 of 1,141,251,880 bytes is 1 part in 42.9. The fraction not held is 0.9767, which is
**42 parts in 43** (0.9767), not 43 in 44 (0.9773).

**Replacement**: "never holds the other 42 parts in 43."

---

## Minor 8 — "what the difference is made of is the allele-length genotype table" is an inference presented as the measurement

**Where**: report: "Take the repeat tracts out — the same run at zero strata — and a sample weighs
**3,391 bytes** instead of 1,141,251. What the difference is made of is the allele-length genotype
table each stratum's fit carries: about **3.35 kB a stratum**."

The measured difference is the whole `StratumFit` record plus its `BTreeMap` node, not only the
genotype vector. At 78 genotypes and a `GenotypeFrequency` of a `SmallVec<[WholeRepeatOffset; 2]>`
plus an `f64`, the vector is on the order of 2.5 kB of the 3.37 kB a stratum; the rest is the
record's own fields (`starts_tried` alone is a `SmallVec<[SlippageStart; 2]>` held inline) and the
map node. The example takes the genotype count as its third argument, so the attribution is one
re-run away from being measured rather than inferred (`-- 338 1 0`).

**Replacement**: "The difference is a whole `StratumFit` record and its map node, about 3.35 kB a
stratum at both counts measured here, of which the allele-length genotype table is the largest
part. Calling never asks for any of it — the length spectrum it does read comes from the cohort
gather, not from these records. Re-running at `-- 338 1 0` would separate the table from the
record."

---

## Minor 9 — "because the seam cannot be written without them" does not hold for the second quantity

**Where**: report: "Two quantities had to be given a route out of the pre-pass first, because the
seam cannot be written without them."

True for the minted-error totals: `assemble` takes the map and asserts it against the error rates.
Not true for the substitution rates — `SsrSampleParameters::by_stratum` is a public field and every
record already carries `substitution`, which the report itself says two paragraphs later ("What was
missing was only the *shape*"). The seam could have folded the records inline.

**Replacement**: "One quantity had no route out of the pre-pass at all and the seam cannot be
written without it; the second had a route and the wrong shape, and giving it the right one keeps
the join out of every consumer."

---

## Minor 10 — the accessor's justification, "one of them was being wrapped into its checked type at each call site"

**Where**: report, *Where the code went its own way*; commit `d8970944`.

At the branch point `ExpectedAlternativeFrequency::try_new` had exactly one non-test call site
outside the type's own tests — `examples/ng_prior_moments_from_reads.rs:861` — plus a test helper in
`run_parameters.rs`, which is **still** wrapping it by hand at HEAD (`run_parameters.rs:686`). The
new accessor has one caller, `from_prepass`. The asymmetry it removes is real; "at each call site"
overstates it.

**Replacement**: "`JointFit::fitted_alternative_frequency`, beside the `fitted_diversity` that was
already there. The two are the pair the genotype prior's seed is built from, and until now only one
of them left the fit already wrapped — the other was wrapped by whoever read it."

---

## Plan conformance, step by step

| step | done? | as specified? |
| --- | --- | --- |
| A1 | yes | yes — field added to `GenericSampleParameters`, filled from `self.minted_errors()` in `estimate` (`estimate.rs:195`), test uses two read groups at different depths and different per-read claims and compares whole-map plus each library's arithmetic. Only the *justification* is wrong (Blocker 1). |
| B1 | yes, smaller | projection not field; premise checked and **the report is right** — `substitution_rates` (`ssr/mod.rs:1188`) already folds every table and omits a key with no compared bases, `assemble_sample_parameters` (`ssr/mod.rs:2249-2258`) already stores the rate on each record and *panics* for a stratum with none. The plan's B1 test as written ("three strata, one with none; the map holds two entries") is unbuildable through the assembled route for that reason — the report should say that outright rather than only citing the three upstream tests. |
| C1 | yes | yes. The message's mechanism is right: `inbreeding` is `None` exactly when the mode is `Fitted` and `rates` has no diploid entry (`estimate.rs:233-235`); `Supplied` always yields `Some` (`fallback.rs:181-187`). So "re-running the fit changes nothing" holds. |
| D1 | yes | seven arguments against the plan's six. Admitted, and the reason holds: no `StratumFits` is carried by `JointFit` or by `SsrSampleParameters`, and every non-test `StratumFits::over` in the tree is a fixture. Types-before-implementation was kept — `d8970944` lands the signature with `todo!()`. |
| D2 | yes | the six reads, the refusal, the call to `assemble`, plus the ownership check the plan did not ask for (an addition, not a divergence). Fixture is three samples over four libraries with the first sample holding two, as specified. |
| E | yes | see the range section below. |

**Divergences the report does not admit:** none of substance. Two omissions worth a line each: the
B1 test could not be written as the plan specified (above), and the seam gained a check the plan
never asked for (`its_own_read_group`) — which the report describes but does not list among the
places the code went its own way.

---

## Scope

Clean. Nothing outside the seam landed: no run driver, no `call_locus` wiring, no fitting of the
batch split. `SequencingBatches` is taken as an argument. The unfitted inbreeding coefficient
refuses by name, with no default and no borrow of `JointFit::hom_excess` — I checked that
`hom_excess` is never read anywhere in the new code. The one addition outside the plan's list,
`JointFit::fitted_alternative_frequency`, is admitted and is a wrap of an existing method.

Nothing inside the scope was skipped.

---

## The range rule

**One sample.** `a_run_of_one_sample_assembles_and_is_uncontaminated` does go through the seam, and
its calibration assertion (0.001 over 0.002 = 0.5) traces to its own input. Two qualifications the
report does not make:

- the fixture *hard-codes* `ContaminationEstimate::NotIdentified { reason: NoPanel }`. That one
  sample yields `NoPanel` is a property of `fit_contamination_over`
  (`joint/contamination.rs:574-576`, `if count < 2 { return refused(NoPanel) }`) — true, but this
  test cannot see it change;
- "it comes back uncontaminated" is a property of `assemble`, already held by the pre-existing
  `a_run_where_nothing_was_identified_is_uncontaminated`. What is genuinely new here is that
  `from_prepass` accepts lists of length one.

Suggested sentence: "**One sample** is a test: `from_prepass` accepts single-element lists and the
calibration comes back traced to that sample's own two numbers. Contamination comes back absent
rather than a fitted zero — which the fixture states, because that a one-sample cohort has no panel
is `fit_contamination_over`'s rule (`count < 2`) and not this seam's."

**A thousand.** Measured, and the shape a sample is given is sourced rather than chosen — that part
is done well. Two things claimed beyond the measurement are Minor 6 and Minor 8 above. The
conclusion the report draws ("a finding for the run driver's plan") follows from the numbers.

---

## Writing

Beyond the replacements already given:

1. **"the seam" is a self-coined name used as if shared.** It appears in the report's second line
   ("the owner set on 2026-08-27: **the seam only**") and about fifteen times after, and the report
   never says in its own words what a seam is — the paragraph that *does* explain it never uses the
   word. CLAUDE.md names this exact failure (*the probe*, *the walk*). Fix: on first use write "the
   one function between the pre-pass's outputs and what calling reads — `RunParameters::from_prepass`
   — which this report calls the seam."

2. **"minted"** does work in five sentences before it is defined. The report's step title glosses it
   once ("The reads' own claim about their errors") and then switches to "minted total",
   "minted-error means", "minted-error total" without ever tying the two together. Add at first use:
   "the *minted* error — the error probability a read carries in from its base and mapping
   qualities, before any fit."

3. **"the accident that has hidden a join in this project four times"** asserts a count with no
   source. Either name where (a report, a commit) or drop the number: "the accident that has hidden
   a join here before."

4. **"three samples over four libraries with no two numbers alike"** (report) is looser than the code
   it describes. `A_RUNS_LIBRARIES` has 0.002 as both an error rate and a minted mean, and 0.004
   twice likewise; the fixture's actual guarantee — which the code's own comment states correctly —
   is that no two values *of the same quantity* are alike. Use the code's wording: "every quantity
   the seam carries differs between every pair it could be swapped across".

5. **`PROJECT_STATUS.md`: "Library target 4,920 → 4,928 passing."** "Library target" is build jargon
   for `cargo test --lib`. Write "The library's own test suite, 4,920 → 4,928 passing."

6. **Report, the `cargo doc` line**: "one link added here pointed at a module path that does not
   exist". The module `ng::calling` exists; the item `ng::calling::RunParameters` does not, because
   the type is not re-exported there. The commit message gets this right and the report does not —
   use the commit's wording.

7. **`minted_errors` as a field name** does not say what the value is (memory: *names must say what
   the value IS*). It is the sum of the reads' own claimed log error probabilities, per read group.
   The type `MintedReadErrors` predates this branch, so this is a note rather than a request:
   `claimed_error_totals` or `reported_error_totals` would carry the meaning without the doc comment.

---

## One structural gap, for the record (Note, not a defect in this branch)

The seam checks that every key a sample carries **belongs to** that sample; it never checks the
converse, that every library of the run **got** a rate. A sample whose highest-numbered library
produced no observations would leave the union's ids contiguous but one short, and
`checked_read_group_count_of` cannot see that — the run then drops that library and panics later in
`ReadGroupParameters::calibration_of` at the first locus carrying one of its reads, which is exactly
the failure the module header (`run_parameters.rs:40-48`) says it refuses at assembly to avoid.
`from_prepass` is the first caller in a position to check it, since it is the first that knows the
run's read-group table. The plan did not ask for it and the report does not claim it; worth an
explicit line in the next plan.
