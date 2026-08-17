# ng cohort merge — C1: one builder's job over one region

*Implementation report, 2026-08-17. Step C1 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[arch](../../ng/arch/cohort_merge.md) §4 and [spec](../../ng/spec/cohort_merge.md) §6.1, §6.2.*

## 1. Plan

Walk, judge, assemble; return the survivors and the failed spans. Everything the walk and the
assembly need already exists (milestones A and B); what C1 adds is **ownership** — which of
the loci a builder can see are that builder's to build.

## 2. The rule, and why it is the whole step

**A locus belongs to the region its first position falls in, and to no other.** The same locus
is closed by every builder whose observations reach it, so without the rule the parallel
arrangement would build it more than once; with it, a builder may *finish* a locus outside its
own region, which is what keeps a locus whole when a deletion carries it past the end
(spec §6.1).

So `build_region` skips a locus starting before its region — an earlier builder sees that one
whole, including the ground it shares with this region — and stops walking at the first locus
starting past its last base, which is a later builder's.

The three verdicts land in three places: `Build` becomes a `CohortObservation`, `Failed`
contributes its span and nothing else, and `TooQuiet` contributes **nothing at all** — ground
the caller examined and found empty, where a failure is ground it refused, and only one of the
two is counted (spec §4.3, §3.3).

## 3. Changes made

In [`build.rs`](../../../../src/ng/run/cohort_merge/build.rs):

- **`RegionOutcome`** — the observations in genome order and the failed spans, with `Default`,
  so that a region which built nothing still delivers one (spec §6.3).
- **`build_region(region, observations_per_sample, max_cohort_locus_span, min_alt_obs)`**.

**Deviation from the architecture's signature, recorded at the code.** Arch §4 takes
`&ObservationCache`; the cache is milestone D. `build_region` takes the slices the cache's
`with_observations` hands out — one per sample, in coordinate order — which is what lets the
serial driver (C2), the oracle every later milestone reproduces, exist before the cache does.
Nothing else in the contract changes: the caller owns the guarantee that the slices reach far
enough for a locus starting inside the region, exactly as the cache will.

**Not implemented, because there is nothing yet to implement them against:** the arch's
"cursor requests stay monotonic" (there are no cursors here — the caller holds the
observations) and "errors end the region and name the sample and span" (nothing in this path
can fail: the walk and the assembly panic on a producer's broken guarantee and return no
errors).

## 4. What the tests pin

6 tests added (104 in the module):

- the three verdicts through one region — one built, one failed span, and a too-quiet locus in
  neither;
- a locus starting **before** the region is the earlier builder's, even though it covers most
  of this region;
- a locus starting **inside** and reaching ten bases past the end comes out whole, with the
  evidence of a sample that only covers its tail;
- a locus starting after the last base is left alone;
- a region with nothing in it still delivers an outcome;
- a locus on another contig is not this region's;
- **and the property the parallel arrangement will rest on**: one walk over 1–100 against ten
  walks over ten bases each gives the same observations, the same alleles, the same per-sample
  support and the same failed spans — with a locus deliberately opening in one ten-base region
  and ending in the next.

**One of my own expectations was wrong again and the code was right**: the fixture yields
three observations, not the two the test's prose claimed — I had forgotten the SNP at 77 while
writing the sentence. Corrected to name all three by position.

## 5. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 104 passed, 0 failed.
</content>
