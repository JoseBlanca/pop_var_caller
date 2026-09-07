# The hidden-duplication filter — C5: the oracles, three of four

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone C, step C5
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §10
**Branch:** `ng-paralog-filter`

## The answer

**Three of spec §10's four oracles hold, on a cohort that produces real records. The fourth cannot
run on this machine: the six tomato accessions' CRAMs are not here.**

| oracle | state |
|---|---|
| (i) off is byte-identical to the pre-filter run | **not run** — needs the tomato slice's reads |
| (ii) on at an unreachable target, stripped, equals off | **holds**, on the mode-equivalence cohort |
| (iii) drop equals tag minus the tagged lines | **holds**, same cohort |
| (iv) both modes agree with the filter on | **holds**, same cohort |

## Where they run, and why there

`a_varying_cohort_on_disk` — the mode-equivalence harness's fixture — is the only cohort in the
crate that **writes records**: at least three, one at each sample's own substitution and one at the
repeat tract the first sample shortens. The command-line fixture (`a_cohort_on_disk`, two samples
and three reads) writes none, which is what C4's review filed as its Blocker.

**So oracle (iv) is what finally covers C4's Blocker.** Each subcommand carries its own copy of the
~45 lines that choose the sink, take it apart, run passes two and three, and print the report.
Running both modes *with the filter on* is the only thing in the crate that compares those two
copies; a divergence between them was invisible to every other test.

## What these oracles do not show

**This cohort is too small for any sample to get a coverage model**, so nothing is ever scored and
nothing is ever flagged. The three oracles therefore prove the *relations between the files* — that
the filter changes nothing it does not flag, that the two modes agree, that tagging and dropping
differ by exactly the tagged records — and **not** that the filter finds duplications. That is the
D milestone's job, on real reads, and the owner has said it will be tested extensively there.

Oracle (iii) is the weakest for this reason: with nothing flagged it asserts only that asking to
tag rather than drop changes nothing else about the file. The stronger half needs a cohort with a
real duplication in it.

## Spec §10's second oracle was amended

**By the owner, 2026-09-07.** As written it could not pass: strip the two `INFO` keys from a
filter-on run and the result still differs from the filter-off file by four header lines — the
`##paralogFilter=` line and the three declarations. Those are emitted only on a run that filtered,
and they have to be: emitting them always would break §10's *first* oracle, that a filter-off run
is byte-identical to the run before this work. The two as written were incompatible. The oracle now
strips the filter's four header lines as well; what it asserts is that **the filter changes no
record it does not flag**, and a header line is not a record.

## What C5 still owes

**Oracle (i), and it needs data this machine does not have.** `benchmarks/tomato1/crams/` is empty
— the reference is on the read-only mount, the reads are not. The standing baseline the plan cites
(2026-09-06, six accessions over the first two 100 kb intervals, 2,311 records, sha256 `84ad19c2…`)
cannot be reproduced or compared against here. Every remaining measurement in the plan — D1's single
accession, D2's six with pass shares, D3's HG002, D4's counts by record kind — waits on the same
files.
