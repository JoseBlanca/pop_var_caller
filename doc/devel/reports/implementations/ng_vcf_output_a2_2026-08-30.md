# ng VCF module, step A2 — what `DP` and `AD` count, decided by the type

**Date:** 2026-08-30. **Branch:** `ng-vcf-output`. **Plan:**
[`../../ng/impl_plan/vcf_output.md`](../../ng/impl_plan/vcf_output.md) step A2.
**Spec:** [`../../ng/spec/vcf_output.md`](../../ng/spec/vcf_output.md) §7.
**Follows:** [step A1](ng_vcf_output_a1_2026-08-30.md).

---

## What landed

`SampleReadCounts` no longer stores the depth. It stores what the depth is *made of* — the
per-allele counts, and the reads no written allele explains — and `DP` is their sum.

**33 tests, all passing** (`cargo fmt --check`, `cargo clippy --all-targets --all-features -D
warnings`, `cargo test --lib ng::vcf`, all green in the container; the runner reports 36 because
the filter also matches three of production's `var_calling::vcf_writer` tests by substring).

## The change, and why it is not just a rename

A1 took `new(allele_reads, depth)` and asserted `depth >= ΣAD`. That assertion is now gone,
because the state it refused **cannot be spelled**: with the components stored and the total
derived, there is no way to hand the type a depth that disagrees with the counts it has to
cover.

The failure this removes is worth naming precisely. A `DP` below its own `AD` is not a crash and
not a parse error — every VCF reader accepts it, and it means something plausible-looking. It
would arise the ordinary way: the worker sums two totals over subtly different sets of reads.
A1's assertion caught it at construction; A2 makes it unrepresentable, which is the stronger
form and the one the plan's principle about silently-wrong numbers asks for.

One arithmetic refusal remains, and it is a different thing: the two counts together must fit a
`u32`. A per-sample depth at one locus is bounded by real coverage — hundreds of reads — so a
total above four billion is a corrupt count, and it is caught rather than wrapped into a small
depth.

## What each column counts — the composition, pinned

Stated against the merge's own fields
([`SampleSupport`](../../../../src/ng/run/cohort_merge/build.rs)):

| | counts | from |
|---|---|---|
| `AD[i]` | reads whose complete observation matched written allele `i`, summed over the sample's read groups | `supported` |
| `DP − ΣAD` | reads the sample observed that no *written* allele explains | see below |
| `DP` | the sum of those two | derived |

**Two things are in `DP` and in no `AD` slot**, and they are why the difference is worth
publishing at all:

- a read whose complete observation matched an allele **candidate selection dropped** — real
  evidence, against a sequence this record does not carry;
- a **partial** observation, a read that ran out inside the locus and so matched nothing exactly
  (`partials`).

**Two more are in neither, and both would be wrong to count as depth:**

- reads that produced **no observation at all** (`reads_without_observation` — every base masked,
  or an `N`). They were never observed in the sense `DP` claims. **And the counter could not be
  used even if the meaning fitted:** the merge sums it per record, so a read silent at two
  records of one locus is counted twice — its own doc says so.
- reads **removed as evidence** (`reads_removed_as_evidence` — named at some of a sample's
  records inside the locus and not at others, so nothing they showed reaches `supported`). No
  usable observation either.

Both are real losses of depth, both are counted by the merge, and neither belongs in a column
that says what the sample showed here. Recording the exclusions is half of what this step is
for: they are the two quantities a later reader would most plausibly fold in by mistake.

## Tests

Three changed or added, all pinning the new structure rather than the old assertion:

- the depth is derived and cannot contradict the allele counts — asserted both as a value and as
  the sum it is defined to be;
- a sample whose reads all miss the written alleles still has depth, with `AD` all zeroes and
  `DP` not — the dropped-candidate case, which is the signal the difference exists to publish;
- the refusal that replaced A1's: two counts too large for the depth column together.

The test that pinned A1's `depth >= ΣAD` refusal was **deleted rather than retargeted**, because
what it asserted is no longer reachable. That is the honest record: a `should_panic` test for a
state the types have made impossible would be testing nothing.
