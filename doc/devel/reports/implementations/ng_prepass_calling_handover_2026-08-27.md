# The pre-pass → calling handover — what was built

**Branch** `ng-prepass-handover`, cut from `main` at `a6e8472b`. Plan:
[`prepass_calling_handover.md`](../../ng/impl_plan/prepass_calling_handover.md), whose scope the
owner set on 2026-08-27: **the seam only**, and **an unfitted inbreeding coefficient refuses the
run**.

## What it is, in one paragraph

The parameter pre-pass measures what a run's data is like and reports it three ways: one value per
sample from the SNP/indel path, one value per sample from the repeat tracts, and one fit over the
whole cohort at once. Calling reads a single object, `RunParameters`. **Nothing built that object
from what the pre-pass produced** — the constructor was called from 29 places at the branch point
and every one of them was inside its own test module, each handing it values written by hand. (The
plan says ten; counted at `629e84ff` it is 29, and the point is the same: none of them was a run.)
There is now one function between the two, `RunParameters::from_prepass` — **this report calls it
the seam**, since it does nothing but join.

**It gathers rather than fits**: every per-library and per-sample number it hands on is one of the
three pre-pass outputs, unchanged. The one derived value is the genotype prior's seed, which
`seed_from_moments` solves in closed form from the cohort fit's two moments and which falls back to
its own constants where a moment is missing.

Two quantities had to be given a route out of the pre-pass first, and one rule had to be decided.
**One of the two had no route at all and the seam cannot be written without it; the other had a
route and the wrong shape**, and giving it the right one keeps the join out of every consumer.

## The library suite

| | tests |
| --- | --- |
| at the branch point (`629e84ff`) | 4,920 passing, 0 failing, 11 ignored |
| now | **4,937** passing, 0 failing, 11 ignored |

`cargo fmt --all -- --check` and `cargo clippy --all-targets --all-features -- -D warnings` both
exit 0. `cargo doc --no-deps --lib` reports **25 unresolved links and exits 101**. It reported 26
before the last commit: one link added on this branch pointed at a module path that does not exist
(`RunParameters` lives in `ng::calling::run_parameters`, not in `ng::calling`) and was corrected.
**The branch point was not measured directly** — 25 is what remains after removing the one link
this branch added, and none of the 25 sites the compiler names is in code this branch wrote.

## Step by step

### The reads' own claim about their errors now leaves the fit (`1eafa835`)

The calling step scores each read at the error probability the read itself claims — its **minted**
error, the number its base and mapping qualities imply before any fit — multiplied by the ratio
between that claim and the rate the pre-pass fitted for its library. It has had the fitted rate for
a while. **The claim, summed per read group while the loci were counted, went out of scope with the
accumulators**, so nothing that assembled a run's parameters could supply it.
`GenericSampleParameters` now carries it, copied from the tally the rest of that value's numbers
were fitted from.

**What made this worth a fixture rather than a line — and the first draft of this paragraph had it
wrong.** It said a missing total takes the defaulted calibration and the run finishes. It does not:
`RunParameters::assemble` requires a read group to have a fitted rate and a minted total or
neither, so an empty map stops the run at assembly naming the first read group, and a map keyed by
another sample's libraries is refused by the seam. **What survives all of that is a map with the
right keys and the wrong numbers in them** — another read group's total under this one's
identifier, which moves every scale it touches and looks exactly like a correct map. So *the map is
not empty* is not what the test asserts. The fixture gives the two libraries different depths and
different per-read qualities, 8 reads a site at 7 nats against 12 at 9, and checks both the whole
map against the tally's own and each library's two numbers against the fixture's arithmetic.

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
fail the new test; **the permutation left every other test of the repeat-tract module green — 201
of them ran**, the 90 in this file among them. (An earlier draft said "in that file", which was the
wrong subject: the file holds 91 tests and the module tree 203, one of them ignored.)

**One thing the plan asked for could not be built.** Its test was to be "three strata, two with
different rates and one with none; the map holds two entries" — and a stratum with no rate cannot
reach a sample's parameters at all, because `assemble_sample_parameters` panics rather than
building a record for it. The rule the plan wanted that test to hold is held three tests upstream,
on the gather itself.

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

The fixture is three samples over four libraries in which **every quantity the seam carries differs
between every pair it could be swapped across** — four error rates, four minted-error means, four
contamination fractions, four tract substitution rates, three inbreeding coefficients, and four
*calibration scales*, which two libraries can share even when their rates and their minted means
both differ. The first sample holds two of the four libraries, so a library's index is not its
sample's; a fixture giving each sample one library makes the two axes the same list of numbers,
which is the accident that hid **five of the seven deliberate defects** the reviews of the calling
loop's contamination step planted, all of which its 4,776 tests passed (PROJECT_STATUS, the
2026-08-26 entry for that step).

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

**The last of those takes explaining, because the naive version does not compile.** The two moments
are different types, so swapping the two arguments is a type error; the mutation had to unwrap each
and re-wrap it as the other, which is what a caller reading them out of the fit by hand could
plausibly do. That the type stops the direct swap is a property worth having and is why
`JointFit::fitted_alternative_frequency` was added rather than the frequency being wrapped at this
call site.

**Two more refusals were added after review**, both for failures that end with a run finishing:
a repeat-tract rate fitted at a ploidy other than the run's, which the lookup can never find, so
every tract would be called on the model's stated constant instead; and a library the run declared
that no sample carries a rate for, which shortens the read-group axis silently and defers the
failure to a locus. The second is reachable from data — a library whose reads were all refused at
admission has no entry anywhere — and **what such a run should get is a design question this does
not settle**: refuse it, as here, or give that library the defaulted calibration a fitted-and-
unmeasurable one would get.

### One sample and a thousand (`2f35e8a4`)

The project requires an answer at both ends of the cohort range.

**One sample** is a test. What is new in it is that the seam accepts single-element lists and the
calibration comes back traced to that sample's own two numbers; the contamination comes back
**absent rather than a fitted zero**, which the fixture states rather than derives — that a
one-sample cohort has no panel is `fit_contamination_over`'s own rule (`count < 2`), not this
seam's. The read likelihood then computes its plain formula, which is the simple case for that
model rather than the weak one. The seam has nothing to special-case.

**A thousand** is a measurement, in `examples/ng_prepass_handover_footprint.rs`, live bytes from
dhat's allocator. It builds the union the way the seam builds it and calls `assemble` directly,
because the only constructor that mints a read-group table without an alignment header is test-only
and an example cannot reach it. At one library a sample holding 338 repeat-tract strata:

| samples | the per-sample results | the run-wide maps | what assembling adds | peak | peak a sample |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1,114,632 | 27,048 | 24 | 1,141,704 | 1,141,704 |
| 10 | 11,146,320 | 265,224 | 240 | 11,411,784 | 1,141,178 |
| 100 | 111,463,200 | 2,659,264 | 2,400 | 114,124,864 | 1,141,248 |
| 1,000 | 1,114,632,000 | 26,595,880 | 24,000 | 1,141,251,880 | 1,141,251 |

**1.14 GB at a thousand samples, and the per-sample results are 97.7% of it** — 1,114,632 bytes
against 26,596 of run-wide maps and 24 of dense vectors. Take the repeat tracts out — the same run
at zero strata — and a sample weighs **3,391 bytes** instead of 1,141,251, so **the repeat-tract
records are 99.7% of what a sample weighs**. That is 3,366 bytes a stratum at 338 strata and 3,334
at 141; running the example with its genotype count set to zero splits it further — **2,496 bytes
of allele-length genotype table** (78 entries of 32 bytes) and 792 bytes of the rest of the record,
plus 79 of the projected rate in the run-wide map. Calling reads one number off each record, the
substitution rate, and nothing else.

**So this is a finding for the run driver's plan and not for the seam.** What a run must hold for
the whole of calling is the 26.6 kB a sample of run-wide maps — 26.6 MB at a thousand. A driver
that projects each sample's substitution rates as that sample finishes, and releases the rest,
peaks at **1/43rd** of the figure above: 26,620 bytes a sample against 1,141,251.

**That 26.6 kB is a lower bound, not the figure.** The measurement leaves contamination empty, so
a run where every library has a fitted fraction carries one estimate per library in the map and one
dense view per library in the result that this table does not count — which is also why its last
column is 24 bytes a library rather than the "vectors" its heading implies. Both are small beside
1.1 MB a sample, and neither changes the conclusion.

The shape a sample is given is stated and sourced rather than chosen: one library, which is every
sample of both benchmark cohorts here; 338 strata a library, from the repeat-tract fit's own report
([`ng_parameter_prepass_ssr_e5_2026-08-13.md`](ng_parameter_prepass_ssr_e5_2026-08-13.md)), against
141 in the tomato SL4.00 catalogue
([`census_tract_grain_b4_2026-08-14.md`](census_tract_grain_b4_2026-08-14.md)) — which measures at
473,484 bytes a sample, the same 3.3 kB a stratum. All three counts are command-line arguments so a
reader can re-run rather than rescale by hand.

**What the measurement does not say.** Nothing in it runs a pre-pass or reads a file; the values
are synthesised at that shape, so it is a footprint measurement and not an accuracy one. And the 78
genotypes a stratum is not a floor over the whole range, as an earlier draft said: the allele
support is `min(repeats, 6) + 7` lengths, so a three-copy tract has 55 genotypes and a four-copy
one 66. It is the value from five copies up, where the support stops growing.

## Where the code went its own way

- **The repeat-tract step is a projection, not a new field**, and the plan's reason for the field
  was checked and is false. See that step above.
- **The seam takes the slippage gather as a seventh argument.** The plan's argument list for it
  names six, and the table one section earlier lists the gather as a pass-through; neither the
  cohort fit nor a sample's parameters carries one, so it has to be passed in.
- **The cohort fit gained one accessor.** `JointFit::fitted_alternative_frequency`, beside the
  `fitted_diversity` that was already there. The two are the pair the genotype prior's seed is built
  from, and until now only one of them left the fit already wrapped in its checked type — the other
  was wrapped by whoever read it. That asymmetry is also what makes the two moments hard to swap by
  accident: they are different types, so the direct swap does not compile.
- **A test-only constructor on `SsrSampleParameters`.** `of_substitution_rates` builds records
  carrying the stated rates and nothing else worth reading, so that the seam's own test does not
  assemble ten fields of a record whose other nine are not the point. It lives beside the type
  whose fields it fills, which is where `ReadGroups::of_libraries` lives for the same reason.
- **The seam checks two things the plan did not ask for.** That every value a sample carries is
  keyed by one of that sample's own libraries, and — after review — that every library the run
  declared got a rate and that no repeat-tract rate was fitted at another ploidy. All three are
  checks rather than values, and each replaces a run that finishes with a message about the run.

## What the reviews found

Three reviews, one brief each: the arithmetic, the tests and their mutations, and plan conformance
with claim-checking. **They found no defect in the seam's arithmetic** — every one of the seven
inputs lands where the doc says, the seed's two moments are in the right order, and no value is
dropped or overwritten on any path a caller can reach. They found **two silent failures the seam
did not check for**, both now refused and both listed above; **eight of 24 planted defects that
every test passed**; and **eleven wrong claims**, of which nine were in this report and its commit
messages and two were doc comments. Four of the eleven were the stated *reason* for a design or a
test, which is the failure this project's review history puts at about 60 in 300.

### The eight tests that could not fail

Six are now fixed and re-checked by re-running the mutation that found them. One is benign and one
is not a defect at all:

| the defect that survived | what it means | now |
| --- | --- | --- |
| each of three ownership checks deleted | the check is called four times and only the error-rate call had a test | four tests, one per route, each naming its quantity |
| the seed's fitted diversity halved | the tracing test read only the seed's mean frequency, which is `α_alt/(α_ref+α_alt)` — exactly `f` for *any* total, so the diversity cancels | the total is asserted against the identity the two moments imply |
| every non-diploid key filtered out of the tract projection | every fixture here is diploid, and the walked one gives all its strata to one library, so two of the key's three axes never varied | four keys differing on each axis in turn |
| two samples sharing a batch | the assertions on those two read the same number, so the axis could not see them exchanged | a batch per sample, and the read-group view checked beside it |
| `fitted_alternative_frequency` returning the heterozygosity | the accessor added here had no test of its own | two, in the module that owns it |
| the ownership check's quantity label crossed | every route's refusal message shares the phrase the test expected | the four expectations name the four quantities |
| `insert` → `entry().or_insert()` | **benign**: the two differ only on a duplicate key, which the ownership check refuses first | left alone, recorded here |
| the repeat-tract length assert pointed at the SNP/indel list | **not a defect**: the assertion above it has already established the two are the same number, so the mutant is equivalent | the new test stays — nothing exercised that assertion before |

**The review's account of that last one is wrong**, and it is worth saying which half: it holds that
a short repeat-tract list would then index past the end of the list and panic with `index out of
bounds`. It would not — the assertion still fires, with the message written for it.

The wrong claims, all corrected in place: that a missing minted-error total lets a run finish;
that a permuted per-sample list silently mis-assigns coefficients and contamination fractions; that
the seam computes nothing; the constructor's caller count (29, not ten); the count and subject of
the tests a mutation left green; two footprint ratios on different bases; the attribution of a
stratum's whole 3.35 kB to its genotype table; "78 genotypes is a floor"; the accessor's symmetry
with its neighbour; and the repeat-tract gather's "refuses to build a record", which is a panic
that aborts the whole sample.

**Four of them were found and fixed before the reviews reported**, by re-checking the report's own
claims against the code: the permuted-list mechanism, the caller count, the doc-link provenance and
the two ratios. The reviews confirmed all four and corrected the caller count again — it is 29, and
the first correction said 28.

## What this does not do, and who owns it

Straight from the plan's out-of-scope list:

- **The run driver** — the program that walks every sample, runs both halves of the pre-pass and
  hands the result to this seam. It is what proves the seam on real data, and the footprint above
  is a finding for it.
- **Wiring `call_locus`** into the merge's builder.
- **Fitting the sequencing-batch split.** The seam takes the batching it is given; the honest input
  until that fit exists is one batch holding everything.
