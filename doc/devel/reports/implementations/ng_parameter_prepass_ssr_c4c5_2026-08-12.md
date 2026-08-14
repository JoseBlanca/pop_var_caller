# ng step 4, the STR path — C4+C5: how a read group's bases read, and the accumulator that files them

*Implementation report, 2026-08-12. Steps C4 and C5 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), taken in one loop, with
the review that followed and the fixes applied — two agents, 28 mutations, 13 behaviour-changing
survivors. Design authority: [`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md)
§2.3, §4 and [`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1.*

**Two plan steps in one loop, named rather than silent.** C4 is the composition channel and C5 is
the front door that uses it; neither is marked *own commit, do not bundle*, and C5's tests are what
make C4's routing checkable at all.

## What the step is

`base_comparison_of` compares each read's bases against the motif **tiled to the length that read
itself shows**, so a read that lost a whole copy is compared against a tract one copy shorter and
its mismatches are substitutions rather than the slip. The two channels then factorise exactly,
which is why the substitution rate is a division and not an axis of any search.

`SsrAccumulators` is the unit's front door: one `StratumTable` per `(read group, stratum)`, one
accumulator per region shard. It borrows each locus and passes it on untouched, and it tallies
rather than repairs — a locus that is not one repeat tract is passed over in silence, and a tract
whose reference length is not a whole number of copies is counted and skipped, never rounded.

**The property C3 could not test is now testable and tested.** A walk cut into shards and merged
equals the uncut walk, entry for entry and counter for counter, over a fixture whose loci are all
deeper than the read cap — below it nothing is drawn and any implementation passes. The reviewer
tried the mutation the plan names, in two shapes (seeding the cap's draw from the first locus a
shard saw, and from a per-shard counter): both die on that test and on nothing else.

## Recorded deviations from the architecture

1. **`base_comparison_of`, not `composition_of`.** It returns a `BaseComparison`, and *base
   composition* is standing genetics vocabulary for nucleotide proportions — the wrong quantity
   for this reader.
2. **`strata()`, `table_for()` and `stratum_count()` rather than a `tables()` map.** Handing out
   `&BTreeMap` publishes the storage decision that `StratumTable`'s own accessor exists to keep
   private.
3. **`new(ploidy)` takes no read-group list**, where arch §4 sketches one. The read groups come
   from each locus's own observations, so a group cannot be missing from a list and silently
   dropped.

## What the review changed

**Blocker — `merge` could drop a whole stratum and nothing noticed.** The branch that adopts a
`(read group, stratum)` key the receiving accumulator has never seen was never reached: every shard
in the fixture held every key at all three cuts. Replacing it with a no-op left the suite green.
In production that deletes exactly the rarest strata — the one a single shard holds, nearest
`MIN_LOCI_TO_FIT`. The fixture now carries a stratum only one shard sees, and a direct test stands
beside it; I reproduced the mutant and it fails both.

**Blocker — nothing checked that a table gets its own library's base counts.** Feeding every
table library 0's comparison, or a zeroed one, both survived: the test that looked like it covered
this compares locus *shapes*, and a shape carries no base counts. Under the first mutant every
library's substitution rate becomes the first library's; under the second every stratum in the
genome reports no rate at all. Neither panics.

**A counter that lied about what it counts.** `loci_subsampled_to_cap` incremented inside the
per-read-group loop, so on a two-library locus it counted two — a field named `loci_…` that can
exceed the loci walked, and on this step's own nine-locus sharding fixture it read 14. The cap does
fire per library, but the counter is named after loci and now counts them.

**Major — the headline claim of the composition channel had no fixture that could fail.** Tiling to
the *reference* tract's length instead of the read's, and comparing against the reference bases
instead of the tiled motif, both survived: for a read shorter than the tract the two rules agree,
the suite's one expanded read tiled perfectly, and every reference tract in the fixtures was pure.
The second mutant is the worse one — it would silently delete the interruption signal the doc's own
caveat says this function charges to the substitution rate.

**A mechanism in my own doc that was wrong, and is now stated with its size.** "Its mismatches are
substitutions rather than the slip" holds for a read whose length differs by *whole copies*. A read
carrying a one-base indel inside the tract — the population the guard bucket holds — tiles out of
phase from the indel onward and is charged about half the bases after it: 11 of 19 in the new
fixture, where the same base lost at the tract's edge costs nothing. Those reads are 9 in 1,000 of
the length-differing reads at a clean stratum and a third to a half of them in the two bands
`GUARD_SHARE_LIMIT` flags, so a stratum above that limit has a substitution rate that should not be
read either. Documented, pinned by test, and worth the owner's eye.

**Also applied:** `// PANIC-FREE:` comments on both narrowings in `base_comparison_of`, with the
saturating `unwrap_or(u32::MAX)` replaced by an `expect` — saturating silently is the one thing
the sibling function in that file argues against; a merge-guard message that no longer claims the
entries were built for different genotype sets (they are not, and the field's own doc says so); a
test for that guard's panic; and a test that a library which witnessed nothing does not cost the
others their entry.

**No wrong numbers this time** — the first step of the six where the numbers check found none.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib parameter_estimation::ssr` | **120 passed** (103 before this step) |
| `cargo test --lib --bins --tests --all-features` | **3,499 passed, 0 failed, 10 ignored** |

Counted rather than recalled: `grep -c '#\[test\]'` gives **29** in `ssr/mod.rs` and **40** in
`locus_offsets.rs`, so the step adds seventeen; the suite moved 3,482 → 3,499. I re-ran three of
the reviewer's mutants against the fixed suite — the dropped stratum, the mixed-up library
comparison, and the per-library cap counter — and each now fails.

**Two gates are red on this branch and neither is this step's**: `cargo clippy --all-targets` fails
in four `examples/` files, and `cargo doc` reports 13 unresolved intra-doc links.

## Audit trail

`tmp/review_2026-08-12_ng-prepass-ssr-c4c5/` — two per-category files and the reviewed patch.
