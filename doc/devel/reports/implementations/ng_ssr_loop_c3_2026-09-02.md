# C3 — the run report partitions what became of the repeat tracts

**Date:** 2026-09-02. **Plan:** [`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md)
Milestone C step C3. **Design:** [`spec/calling_loop_ssr.md`](../../ng/spec/calling_loop_ssr.md)
§3.2. **Modules:** [`src/ng/run/callers.rs`](../../../../src/ng/run/callers.rs) (`TractOutcomes`)
and [`src/ng/run/report.rs`](../../../../src/ng/run/report.rs).

---

## What landed

**One count became five, and the five are a partition.** `tract_loci_set_aside` — a single number
that meant *built and not scored* — is replaced by `TractOutcomes`: called, refused as
`notPeriodic`, called over fewer sequences than segregate (`tooManyAlleles`), refused for a
candidate carrying no whole motif copy, and set aside as a repeat cluster. Every tract-kind locus
the merge built is in exactly one of them, and `built()` is their sum.

On a real run — HG002 at 30× over 200 Tier intervals — the report now opens the tract section
with:

```
repeat tracts: 24 built, of which 24 called
```

and adds a line per refusal only where one happened.

## Why the refusals have to be counted: the file cannot say them

**A tract refused as `notPeriodic` leaves no record.** Selection narrows it to the reference tract
alone, so every sample is called homozygous reference, the locus establishes no variant, and it is
left out of the file — where it is indistinguishable from a tract nobody varied at
([`vcf_output.md`](../../ng/spec/vcf_output.md) §9). The same holds for the two that are never
scored at all. So three of the five outcomes are visible **only** in this count, which is what
makes it worth a line rather than a footnote.

**A truncated tract is counted apart from the cleanly called ones**, although it *is* called —
over the sequences the cap kept. Folding the two together would give a reader one number they
would take for clean calls.

## Each line appears only when it happened

Four zeroes under a headline is a report a reader stops reading, so a refusal prints only where
its count is non-zero — and the whole section prints only when the run built a tract at all. Both
are asserted, in both directions.

## Tests — 3 new, 1 rewritten

| test | what it pins |
|---|---|
| `tract_loci_the_run_could_not_score_are_a_line_of_their_own_and_only_when_there_are_some` (rewritten) | the headline's sum, a line per refusal with its own wording, and that a run where a refusal never happened prints no zero for it |
| `the_tract_outcomes_sum_to_the_tracts_built` | the five are a partition, and what `refused_by_a_filter` adds up |
| `each_selection_verdict_lands_in_its_own_tract_outcome` | the verdict-to-outcome mapping, including that a truncated tract is **not** counted among the cleanly called |

**The third test exists because a mutation predicted to survive did.** Counting a `notPeriodic`
tract as called changed nothing anywhere: no fixture in the suite builds a non-periodic tract
through the driver, and building one end to end is a fixture this step does not need. Extracting
the mapping into `count_this_tract` made it reachable, and the test kills the mutation.

## Mutation testing — three run, three killed

| mutation | outcome |
|---|---|
| the headline prints even when the run built no tract | killed |
| `built()` forgets the bundles | killed |
| a `notPeriodic` tract is counted as called | killed — **by the test added after it survived** |

## Validation

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` 6,013 passed / 0 failed / 14 ignored, all in the container. The report above is
from a real run of the shipped binary, not a fixture.
