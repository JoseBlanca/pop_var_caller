# Where `estimate-parameters` memory goes: more on its progress lines

*Implementation, review and fixes, 2026-09-26. Branch `estimation-instrumentation`. Asked for by
the owner after the simulator could not reproduce kimura's 100-sample figures
([plan](../../implementation_plans/estimation_memory.md) §7): instead of calibrating a simulation,
run the merged build on kimura and have it say what each step holds.*

## What the progress lines now add

Each line still ends in the process's resident memory and the most it has held (`VmRSS`, `VmHWM`).

| line | new information |
|---|---|
| `censuses and reference read` | how many samples, and how many read groups they hold |
| `SNP/indel fit: evidence read` | the bytes the ordinary-position evidence holds, in total and per read group; how many non-reference observations; the samples' mean depth (lowest, median, highest, from the stored codes after the cap); how many samples have no walked position |
| `contamination:` (new) | how many markers of how many positions, over how many units a fraction is fitted for, and the bytes the markers hold |
| `repeat-tract evidence: arranged` (new) | tracts the guard left out; (tract, sample) pairs with reads out of the pairs the kept tracts make, as a rate in 100; the bytes the arranged evidence holds; bytes a pair with reads |

The byte figures are the vectors' reserved sizes (`GenericEvidence::heap_bytes`,
`StratumEvidence::heap_bytes`, the markers' vectors); the allocator rounds each block up, so the
process holds somewhat more. The difference between these and the resident memory is itself
something the run will show.

These answer the two questions the simulator could not: how many read groups a real sample holds
and what one costs, and at how many tracts a real sample has reads.

## What it must not change

Nothing but stderr. No computed value is fed back; the depth summary sorts a copy of the
per-sample depths. The cross-platform checksum test passes unchanged.

## Tests

`progress::tests::a_size_reads_in_mebibytes_below_a_gibibyte`,
`census::tests::a_read_groups_evidence_counts_its_codes_and_its_observations`,
`ssr_fit::the_largest_table::a_stratums_rows_and_bytes_are_counted_over_every_tract`.

## Review

One reviewer in its own worktree, all applicable categories: no Blocker, no Major; nothing can
change a result or panic (empty cohort, zero read groups and NaN depths checked), and the added
walks cost seconds at most at 2,000 samples. Applied: the pair count now says it leaves out the
tracts the guard dropped (Minor); samples with no walked position are counted apart, since their
depth stands in at 1.0 (Minor); contamination's count is of *units*, a library each or a sample
each at the sample grain (Nit); the markers' own struct size is counted (Nit). Not applied: two
style Nits (a one-line forwarder `table_size`, and long paths to `progress::size`).

## Validation

In the dev container: fmt and clippy `-D warnings` clean; targeted tests (census_fit,
contamination, estimate_parameters, cross_platform and the three new ones) 73 passed; the full
suite before the review's fixes 4,879 passed, 3 failed — the pre-existing failures in
`examples/ng_generic_loci_dump.rs` and `examples/ng_ssr_loci_dump.rs`.
