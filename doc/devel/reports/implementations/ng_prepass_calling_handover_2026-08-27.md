# The pre-pass → calling handover — what was built

**Branch** `ng-prepass-handover`, cut from `main` at `a6e8472b`. Plan:
[`prepass_calling_handover.md`](../../ng/impl_plan/prepass_calling_handover.md), whose scope the
owner set on 2026-08-27: **the seam only**, and **an unfitted inbreeding coefficient refuses the
run**.

## What it is, in one paragraph

The parameter pre-pass measures what a run's data is like and reports it three ways: one value per
sample from the SNP/indel path, one value per sample from the repeat tracts, and one fit over the
whole cohort at once. Calling reads a single object, `RunParameters`. **Nothing built that object
from what the pre-pass produced** — the constructor was called from 28 places at the branch point
and every one of them was inside its own test module, each handing it values written by hand. (The
plan says ten; counted at `629e84ff` it is 28, and the point is the same: none of them was a run.)
There is now one function that does it,
`RunParameters::from_prepass`, and it computes nothing: every number it hands on was measured by
one of the three.

Two quantities had to be given a route out of the pre-pass first, because the seam cannot be
written without them, and one rule had to be decided.

## The library suite

| | tests |
| --- | --- |
| at the branch point (`629e84ff`) | 4,920 passing, 0 failing, 11 ignored |
| now | **4,928** passing, 0 failing, 11 ignored |

`cargo fmt --all -- --check` and `cargo clippy --all-targets --all-features -- -D warnings` both
exit 0. `cargo doc --no-deps --lib` reports **25 unresolved links and exits 101, exactly as it did
at the branch point** — one link added here pointed at a module path that does not exist and was
corrected, which is what took it back to 25.

## Step by step

### The reads' own claim about their errors now leaves the fit (`1eafa835`)

The calling step scores each read at the error probability the read itself claims — from its base
and mapping qualities — multiplied by the ratio between that claim and the rate the pre-pass fitted
for its library. It has had the fitted rate for a while. **The claim, summed per read group while
the loci were counted, went out of scope with the accumulators**, so nothing that assembled a run's
parameters could supply it. `GenericSampleParameters` now carries it, copied from the tally the
rest of that value's numbers were fitted from.

**What made this worth a fixture rather than a line.** A read group with a fitted rate and no
minted total is not refused downstream. It takes the defaulted calibration — scale one, every read
of the library charged the error floor — and the run finishes normally with quietly overconfident
reads. So *the map is not empty* is not what the test asserts. The fixture gives the two libraries
different depths and different per-read qualities, 8 reads a site at 7 nats against 12 at 9, and
checks both the whole map against the tally's own and each library's two numbers against the
fixture's arithmetic.

Emptying the map and swapping the two libraries' totals both fail it. Both left the other ten tests
in that file green — including the one asserting that the two entry points return equal values,
because they share the code that was broken.

### A sample's tract substitution rates, in the shape calling takes (`ebcd86a8`)

A read out of a repeat tract can disagree with the allele it came from in two ways: it can report
the wrong number of repeats, which is slippage, or the right number with a wrong base inside, which
is the substitution rate. Calling scores those separately and needs both keyed by stratum.

**⚑ The plan's premise for this step turned out to be false, and the step is smaller than it
planned.** The plan said the gather did not exist and that `SsrSampleParameters` "carries the fits
and drops the tables, so the two halves never meet". Both halves are in fact already there:
`substitution_rates` folds every stratum's table into exactly that map and already omits a stratum
that compared no bases rather than calling it zero — the rule the plan wanted, with three tests
holding it — and `assemble_sample_parameters` stores each rate on the record it builds, unchanged.
What was missing was only the *shape*: a map keyed by stratum, which a consumer would otherwise
have to join for itself.

So what shipped is `SsrSampleParameters::substitution_rate_by_stratum`, a **projection** of the
records rather than a second field beside them. A field would have been a second copy of numbers
already on the records, free to drift from them.

The fixture gives its three strata unlike rates — 2, 5 and 10 mismatched bases in a thousand —
because every other fixture that reaches these records is built from reads that show their tract
perfectly, so every rate is the same measured zero and a rate read from under a neighbour's key is
the right number. Pairing each key with its neighbour's rate, and keeping only the first key, both
fail the new test; **the permutation left all 201 other tests in that file green.**

### A sample with no inbreeding coefficient stops the run, by name (`5f61bed4`)

`ParameterEstimationError::InbreedingNotFittedForSample`, beside the three inbreeding failures
already in that enum, and it has to be tellable apart from them. Those three are a search that ran
and did not settle, and the answer is to widen it or supply the number. This one is a search that
was never run: the coefficient is measured on the **diploid** part of a genome, and a sample with
no diploid region reports none without that being a failure — above two copies the quantity needs
several identity-by-descent coefficients and is deferred, below two there are no heterozygotes to
be short of. So re-running the fit changes nothing, and the message says so, so that a reader does
not go and re-run it.

**No default, and in particular not the cohort fit's homozygote excess** — which is available and
is the obvious thing to reach for. It is measured by the very fit whose diversity the coefficient
exists to correct, so using it makes the correction circular; the census-moments report already
carries that warning. The rule the pre-pass states one level up is the same one: a cohort's
diversity divides by `1 − F`, so a coefficient invented rather than measured is amplified rather
than absorbed.

### The seam (`d8970944` types, `9ac41520` implementation)

`RunParameters::from_prepass` takes the two per-sample lists, the cohort fit, the run's read-group
table, the declared sequencing batching, the repeat-tract slippage gather and the ploidy, and
returns a `RunParameters` or the refusal above.

**Two joins that no type enforces, and both are silent when wrong.**

- **Sample order.** The two per-sample lists are joined to the run's sample table **by position**.
  Sample order is therefore an argument rather than a map's iteration order, which is the one
  design decision in this step.

  **The inbreeding coefficient is the one value with no identifier on it.** Every other quantity a
  sample carries is keyed by read group, and the cohort fit's per-sample results are keyed by
  sample name — so those either land under a key that says whose they are, or are looked up by the
  name the read-group table gives that position. A coefficient is a bare number, and a permuted
  list sends it to the wrong sample with nothing about the value saying so.
- **The read-group union.** The run's per-library maps are the samples' maps put together, which is
  a union only because a read group belongs to exactly one sample: the run's read-group table files
  each declared `@RG` under the single sample its header names. So every key a sample carries is
  checked against that sample's own read groups. A value that fails it would land under a real
  identifier belonging to somebody else, and every read of that library would then be scored under
  another sample's chemistry, with nothing downstream able to tell.

  **That check is also what catches a permuted list**, which the first draft of this report and of
  the seam's own doc comment got wrong: they said a permuted list "gives every sample a coefficient
  and a contamination fraction, just not its own", and neither half holds. Contamination is looked
  up by the name at that position and cannot be permuted at all. And a sample handed its
  neighbour's results carries its neighbour's read-group identifiers, so the check fires and the
  run stops rather than proceeding with a wrong coefficient. The residual gap is one sample
  carrying *no* read-group-keyed value at all, which the SNP/indel fit does not produce —
  `GenericAccumulators::estimate` refuses a sample with no read group with reads.

The fixture is three samples over four libraries with **no two numbers alike** — four error rates,
four minted-error means, four contamination fractions, four tract substitution rates, three
inbreeding coefficients, and four *calibration scales*, which two libraries can share even when
their rates and their minted means both differ. The first sample holds two of the four libraries,
so a library's index is not its sample's; a fixture giving each sample one library makes the two
axes the same list of numbers, which is the accident that has hidden a join in this project four
times.

**The four substitution rates share one stratum on purpose.** With a stratum apiece, a rate swapped
between two libraries lands on a key nothing asks about and the lookup answers *absent*, which is a
visible failure. Sharing the stratum makes a swap answer with another library's number.

Five deliberate defects, all killed:

| the defect | what caught it |
| --- | --- |
| the samples read in reverse order | the tracing test, and the refusal test |
| every sample's contamination looked up under the first sample's name | both, again |
| every sample's tract rates taken from the first sample's | the tracing test |
| the library-ownership check disabled | the refusal test written for it |
| the seed's two moments swapped | the tracing test, and the one-sample test |

### One sample and a thousand (`2f35e8a4`)

The project requires an answer at both ends of the cohort range.

**One sample** is a test: the run assembles, and it comes back **uncontaminated** — absent rather
than a fitted zero — because there is no panel for a stray-read fraction to be surprised by. The
read likelihood then computes its plain formula, which is the simple case for that model rather
than the weak one. The seam has nothing to special-case for it.

**A thousand** is a measurement, in `examples/ng_prepass_handover_footprint.rs`, taken with dhat's
allocator as live bytes. At one library a sample holding 338 repeat-tract strata:

| samples | the per-sample results | the run-wide maps | what assembling adds | peak | peak a sample |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1,114,632 | 27,048 | 24 | 1,141,704 | 1,141,704 |
| 10 | 11,146,320 | 265,224 | 240 | 11,411,784 | 1,141,178 |
| 100 | 111,463,200 | 2,659,264 | 2,400 | 114,124,864 | 1,141,248 |
| 1,000 | 1,114,632,000 | 26,595,880 | 24,000 | 1,141,251,880 | 1,141,251 |

**1.14 GB at a thousand samples, and 99.7% of it is nothing the seam reads.** Take the repeat
tracts out — the same run at zero strata — and a sample weighs **3,391 bytes** instead of
1,141,251. What the difference is made of is the allele-length genotype table each stratum's fit
carries: about **3.35 kB a stratum**, at both stratum counts measured here. Calling never asks for
it.

**So this is a finding for the run driver's plan and not for the seam.** What a run must hold for
the whole of calling is the 26.6 kB a sample of run-wide maps — 26.6 MB at a thousand. A driver
that projects each sample's substitution rates as that sample finishes, and releases the rest,
never holds the other 43 parts in 44.

The shape a sample is given is stated and sourced rather than chosen: one library, which is every
sample of both benchmark cohorts here; 338 strata a library, from the repeat-tract fit's own report
([`ng_parameter_prepass_ssr_e5_2026-08-13.md`](ng_parameter_prepass_ssr_e5_2026-08-13.md)), against
141 in the tomato SL4.00 catalogue
([`census_tract_grain_b4_2026-08-14.md`](census_tract_grain_b4_2026-08-14.md)) — which measures at
473,484 bytes a sample, the same 3.3 kB a stratum. All three counts are command-line arguments so a
reader can re-run rather than rescale by hand.

**What the measurement does not say.** Nothing in it runs a pre-pass or reads a file; the values
are synthesised at that shape, so it is a footprint measurement and not an accuracy one. And the 78
genotypes a stratum is a floor — that is a five-copy dinucleotide tract at two genome copies, and
longer tracts carry more.

## Where the code went its own way

- **The repeat-tract step is a projection, not a new field**, and the plan's reason for the field
  was checked and is false. See that step above.
- **The seam takes the slippage gather as a seventh argument.** The plan's argument list for it
  names six, and the table one section earlier lists the gather as a pass-through; neither the
  cohort fit nor a sample's parameters carries one, so it has to be passed in.
- **The cohort fit gained one accessor.** `JointFit::fitted_alternative_frequency`, beside the
  `fitted_diversity` that was already there. The two are the pair the genotype prior's seed is built
  from, and one of them was being wrapped into its checked type at each call site — which is how
  two numbers that travel together come to be wrapped differently.
- **A test-only constructor on `SsrSampleParameters`.** `of_substitution_rates` builds records
  carrying the stated rates and nothing else worth reading, so that the seam's own test does not
  assemble ten fields of a record whose other nine are not the point. It lives beside the type
  whose fields it fills, which is where `ReadGroups::of_libraries` lives for the same reason.

## What this does not do, and who owns it

Straight from the plan's out-of-scope list:

- **The run driver** — the program that walks every sample, runs both halves of the pre-pass and
  hands the result to this seam. It is what proves the seam on real data, and the footprint above
  is a finding for it.
- **Wiring `call_locus`** into the merge's builder.
- **Fitting the sequencing-batch split.** The seam takes the batching it is given; the honest input
  until that fit exists is one batch holding everything.
