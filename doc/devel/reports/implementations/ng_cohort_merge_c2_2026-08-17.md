# ng cohort merge — C2: the serial driver, and the oracle

*Implementation report, 2026-08-17. Step C2 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §15, and Checkpoint C's own claim.*

## 1. Plan

One builder, no cache, no organiser, no threads, over the whole analysed stretch. **Its
output is what every later milestone must reproduce**, which is the only assertion that can
catch an ownership or overlap defect once the work is divided (spec §15).

## 2. What it is, and why it is thin

[`merge_cohort_serially`](../../../../src/ng/run/cohort_merge/serial.rs) walks the run's
analysed regions in order, calls [`build_region`] on each with the whole per-sample slices,
and gathers the observations and the failed spans.

With the whole stretch in hand there is nothing else to do: which loci a builder closes
depends on which observations it can see, and here it sees all of them. What the later
milestones add is a **narrower view** per builder — a window from the cache (D), regions in
parallel (E) — and the loci that differ between the two views are exactly what the organiser's
overlap resolution exists to settle. That is what makes this the oracle: it has no window and
no resolution to get wrong.

**No overlap resolution happens here and none is needed**, for the same reason C1's review
established for `build_region`: every region is built from the same complete view, so the loci
closed over one region are the loci closed over any other, and the ownership rule alone keeps
each exactly once.

## 3. Changes made

- New file [`serial.rs`](../../../../src/ng/run/cohort_merge/serial.rs) — the driver and its
  tests. **A file of its own**, because driving is a different job from assembling: `build.rs`
  answers "what is this locus", and this answers "which ground is walked, in what order".
- One line in `mod.rs` to declare it.

## 4. The end-to-end fixture Checkpoint C asks for

`a_cohort_observation_is_built_from_minted_observations` is the milestone's own claim: **two
samples' reads on disk, through the real generic locus generator, into one cohort
observation.** Two BAMs are written to temporary directories, opened as samples, and walked by
`PileupGenerator` over chr2 of the shared fixture reference (all `A`):

- **sample 0** carries a substitution at position 112 in three 30-base reads;
- **sample 1** carries a five-base deletion over 110–114 in three 35-base reads.

They chain into **one** cohort locus at 109–114 whose alleles are `AAAAAA` (the reference over
six bases), `AAACAA` (the substitution widened onto them) and `A` (the deletion).

**What that exercises that no fabricated fixture had:** the generic mint writes a record at
**every covered position**, so inside this six-base locus the substitution sample has **six**
records and each of its three reads is named at all six — its allele is composed across them
and its quality sums are divided, which is the case B3's own review found untested. The test
asserts the composed quality against the six records' own means (the weakest, times three
reads) rather than against a constant, and the deletion sample's support against the record
the generator wrote, so the "one record is exact" and "several records divide" halves are both
pinned on minted data.

**Two things the fixture taught:**

- **The 30-base minimum read length is a real gate.** The deletion reads were 25 bases of
  sequence at first and the filter dropped every one, so the sample minted nothing at all —
  visible only because the test asserted that both samples minted something before merging.
- **A gap between a sample's records means no coverage, not "nothing to say".** B2's
  `alleles_of_sample` doc explains the reference-filling between two of a sample's records as
  ground "where this sample minted nothing, because none of its reads departed from the
  reference there". With a per-position mint that is not what a gap means — and the case
  cannot arise on the generic path at all, because a read named at two records of a
  per-position mint is named at every record between them. The doc is corrected in this
  commit.

## 5. What the other tests pin

- the analysed regions are walked in order and a locus reaching from one into the next is
  built once, by the first;
- several contigs come out in the order they were analysed;
- the failed spans are gathered across regions while a too-quiet locus reaches neither vector;
- analysing nothing yields nothing rather than failing.

## 6. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 111 passed, 0 failed.

## 7. What C2 does not do

- **No cache and no threads** — milestones D and E, and the partition-invariance assertion
  against this oracle is E4's.
- **It is quadratic if it is misused.** A builder pays for every observation preceding its
  region (3.3 µs per base at 63 samples, C1's review), so handing this driver thousands of
  short regions with the whole slices would be quadratic in the stretch. It is written for
  whole analysed regions, where that prefix is empty, and the cache is what makes short
  regions affordable.
</content>
