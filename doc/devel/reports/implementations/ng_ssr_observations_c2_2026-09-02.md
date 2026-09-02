# ng STR observations — C2: the tract slot is filled, and the driver sets each locus aside

*2026-09-02. Step C2 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §3.1 and §5](../../ng/spec/run_ssr_observations.md). Branch `ng-ssr-observations`.*

## Plan

The tract generator was built, tested and reachable only from the development tools; a run's
generator set refused its slot as unbuilt. Filling it makes a run *produce* repeat-tract
observations, per sample and then per cohort. What the calling loop does with one belongs to
the other plan, so the driver sets each aside and counts it.

## Assumptions

**The bundle radius comes from the criteria the ground was cut with**, not from the constant
that happens to equal it. The generator's constructor checks that the flank it fetches fits
inside that radius, and checking it against a constant would make the check answer a question
nobody asked. No flag moves the radius today; if one is added, this is the second place it has
to reach, and the code says so.

**The set-aside is counted, not collected.** A run over tract-rich ground meets millions of
these and what a report owes is how many, not where each was.

**It is a new count rather than an existing one.** `unhandled_not_implemented` counts *regions*
whose generator slot is unfilled — bundles, now that the tract slot is filled. This counts
*loci* that were built, merged across the cohort, and then not called. Merging them would make
the run report unable to say which of the two happened.

## Changes made

| file | change |
|---|---|
| `run/walker.rs` | `generic_path_generators` builds an `SsrGenerator` into the tract slot, with its own reference accessors and the run's bundle radius |
| `run/mod.rs` | `RunError::TractGeneratorSettings`, transparent over the generator's own config refusal |
| `run/callers.rs` | `set_aside_because_nothing_calls_it_yet`; both drivers count and skip; `tract_loci_set_aside` on `CalledCohort` and `WrittenCohort` |

The aligner is not a knob: the delimiter bake-off is recorded and a run gets its winner, the
unit-robust algorithm 4u.

## Tests added and changed

**Added** `a_tract_a_sample_varies_at_is_built_and_set_aside_uncalled` — one sample whose reads
span the tract with the fifteen bases of flank the generator fetches on each side, carrying one
changed base inside it, so the merge has something to build. It asserts **both halves**: the
locus is counted as set aside, *and* no called locus lies over the tract's ground. Asserting
only the count would pass on a run that set the locus aside and called it as well.

**Changed** the E2 thread-invariance test, in two places, and both are the behaviour moving
rather than the test being loosened:

- each sample's `unhandled_not_implemented` over that fixture's three segments was **1** and is
  now **0** — the tract was a region whose generator did not exist and is now a region with one;
- `tract_loci_set_aside` joins the per-thread-count comparison, so the new count cannot depend
  on the pool size either.

**And the mixed cohort's own set-aside count is zero**, which is worth stating because it looks
like the guard not working: no sample of that cohort varies inside the tract, so the merge finds
it too quiet to build and there is no locus to set aside. The assertion says so and names the
fixture that does vary there.

## Validation

In the dev container. `cargo fmt --check` and `cargo clippy --all-targets --all-features -D
warnings` clean. `cargo test --lib --tests --examples --all-features --no-fail-fast` — **5,931
passed, 0 failed, 14 ignored** in the library suite (5,929 before this step); every integration
target green; the three known locus-dump failures and the psp writer bench unchanged.

**Two mutations, both killed:**

| mutation | what a run would do | tests failed |
|---|---|---|
| the guard never fires | score a tract under the SNP/indel model and write a plausible record | 1 |
| the slot goes back to unfilled | build no tract observation at all | 2 |

## Tradeoffs and follow-ups

- **A tract observation is built and then thrown away**, which is work a run pays for and gets
  nothing back from until the calling loop's dispatch lands. That is the price of the two plans
  being independently buildable, and it is bounded: the routing change of B1 cut the ground
  this happens on by about seven times.
- **The bundle slot stays unfilled.** A bundle has neither a generator nor a design.
- **The run report does not print the new count yet** — C3.
