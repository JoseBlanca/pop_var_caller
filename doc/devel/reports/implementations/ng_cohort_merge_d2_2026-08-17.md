# ng cohort merge — D2: the merge, read through the cache

*Implementation report, 2026-08-17. Step D2 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §6.4, with §15's oracle as the thing to reproduce.
Revised after [the D2 review](../reviews/ng_cohort_merge_d2_2026-08-17.md), which corrected
four claims in the first draft — each is fixed in place and named in §7.*

## 1. Plan

Drive the merge through [D1's cache](ng_cohort_merge_d1_2026-08-17.md) instead of holding every
sample's observations at once, and prove the output is C2's.

## 2. What it is

[`merge_cohort_through_cache`](../../../../src/ng/run/cohort_merge/serial.rs), beside the
oracle it must reproduce. Per building region: evict at its first base, cover it, hand the
window to `build_region`, gather. The oracle stays exactly as it was — **it is the thing being
reproduced, so it is not touched**, apart from the guard both drivers share (§5).

Three things came out of writing it, and all three are in the interface:

- **The cache is handed in, not made here.** It holds one reader per sample *for the whole run*
  (spec §6.4), so it outlives one merge, and at milestone E it belongs to the organiser. It is
  also the only way to observe what this driver does to memory (§6).
- **The division into building regions is its own function**,
  [`building_regions_of`](../../../../src/ng/run/cohort_merge/organise.rs) — in the organiser's
  file, because handing regions out is what the organiser does at milestone E and a second
  derivation of the clamp there is exactly what its tests exist to catch.
- **`ObservationCache::held_observations_len`** reports the window's size, for §6.

## 3. The claim, and how it is proved

**Byte-identical to `merge_cohort_serially` on the same observations.** The comparison is the
`Debug` rendering, locus by locus — `CohortObservation` has no `PartialEq`, a field-by-field
comparison would silently stop covering a field added later, and two distinct `f64` sums render
as distinct strings, so a quality divided differently shows. Rendering entry by entry rather
than whole is the review's: on failure the whole-outcome form printed 36 kB with no pointer at
the difference, against 888 bytes naming the first entry that differs.

| fixture | what it puts in the way |
|---|---|
| 60 single-base loci, a 26-base deletion, and a sample inside it, at widths 1, 3, 20, 47 and 600 | building regions narrower than a locus, and one as wide as the whole stretch |
| a record at 55 against an analysed region ending at 50 | the last building region's clamp |
| two analysed regions, two refused loci, at 7-base regions | failed spans and quiet ground crossing several building regions |
| two contigs, a record at 900 on the second | a reader crossing a contig, and eviction at the next contig's first base |
| **two samples' reads on disk, through the real generic generator**, at 4-base regions | Checkpoint C's own locus, at a width narrower than the six bases it spans |
| **200 random layouts** — 1 to 6 samples, 1 or 2 contigs, records up to 60 bases, analysed regions with gaps, widths 1 to 30, several span bounds and keep thresholds | everything the five above do not enumerate |

The randomised differential is the review's, which ran it at 600 layouts with **no
disagreement**; it is kept because it killed two mutations the fixtures did not, and because
every one of its layouts contains a record that opens inside an analysed region and ends past
it — the regime where the cache must draw past the analysed ground and the oracle simply
already has it.

## 4. What the division has to get right

**The building regions tile each analysed region exactly and the last one is clamped to it.** A
locus belongs to the builder whose region holds its first position, so a last region running to
its full width would claim loci past the analysed ground — which the oracle never builds.

At the top of the coordinate space the division stops rather than wrapping: the successor is
`checked_add`ed, and there is none. Its test takes three regions rather than collecting the
iterator, because the failure it guards against is an *unbounded* division — an unbounded
collect would take the test binary down instead of printing a difference.

## 5. The hole the review found in the claim itself

**An analysed region whose ends are the wrong way round was read differently by the two
drivers**, and nothing refused it: the division orders the two ends, `build_region` does not.
On `50-1` over records at 12 and 45, the oracle built **nothing** and the cached driver built
**both loci** — and with a second region after it, the same locus came back **twice**, which is
the corruption the neighbouring guard exists to prevent. The shared check now refuses a region
whose ends are inverted, and the guard is named for the two rules it enforces
(`refuse_malformed_analysed_regions`).

A second gap in the same guard: both refusal fixtures overlapped by 21 bases, so relaxing its
comparison from `<` to `<=` passed every test while duplicating any locus on a shared base. Two
regions sharing exactly one base now pin it.

## 6. Two properties the merge's output cannot show, and what pins them

Both were found by mutating the code rather than by reading it, and the review found the first
pin insufficient.

- **That the ground is divided at all.** A driver that ignored the width and built each analysed
  region as one produces *exactly the right answer* — the oracle does precisely that. Testing
  `building_regions_of` directly is not enough: nothing said the driver *called* it, and
  mutating the driver's own call site left the division test green. What pins it now is
  `the_width_the_caller_asks_for_is_the_width_the_driver_builds_in`, which merges the same
  fixture at two widths and reads the cache: **60 observations held at one region for the
  stretch, 2 at twenty-base regions.**
- **That the driver evicts.** Same shape: never evicting gives the same observations while
  holding the whole stretch. `the_driver_evicts_as_it_goes_and_the_window_stays_short` asks the
  cache afterwards — 2 held against the 60 the stretch holds — and
  `the_window_stays_short_up_to_a_failure` reads it *mid-stretch*, at a source failure, which is
  the only place from outside where the driver's drawing pace is visible.

## 7. Five claims in the first draft of this report were wrong

All were the author's own about the author's own fixtures, and every figure quoted from another
document was right — the pattern the plan's skill names:

- **"14 new tests" was 12** (`git diff | grep -c '#\[test\]'`), and "every one of them can fail"
  was a different claim from its evidence: nine mutations of the driver, each killing at least
  one test, does not establish that each test is killable.
- **"the same deletion at 20-base regions" was listed as asserting byte-identity** — that test
  runs the cached driver alone. The fixture that does compare there is the clamp one.
- **"the window holds 2, being the last region's record and the one draw past it"** — the two
  are the records at **581 and 591, both inside the last building region 581–600**, and the
  source is spent after 591, so nothing was drawn past the analysed ground. On the same fixture
  at five-base regions the window ends holding **none**: the number moves with the width and the
  record spacing rather than being a bound a forward reader pays.
- **"the width ignored kills 2 tests"** — true of a mutation inside `building_regions_of`, and
  false of one at the driver's call site, which is the shape the defect takes. That is finding §6.

## 8. Tests — 22 new across the two files, and the mutations that check them

`serial.rs` goes from 9 tests to 26 and `organise.rs` from 27 to 32 (two of those five moved
with `building_regions_of`); the module's suite is 164.

Thirteen mutations were written against the driver, the division and the accessor, and the
module's tests re-run in the container for each: **all thirteen fail at least one test.** Five
of the thirteen are mutations the reviewers found surviving the first version. The driver is
`tmp/mutate_d2.sh` (scratch, not tracked).

| mutation | killed by |
|---|---|
| the overlap guard accepts one shared base | `analysed_regions_sharing_one_base_are_refused` |
| an inverted analysed region is accepted | `an_analysed_region_with_inverted_ends_is_refused` |
| the last building region is not clamped | 4 tests |
| the width is ignored at the driver's call site | 3 tests |
| the cover covers the whole analysed region | `the_window_stays_short_up_to_a_failure` |
| eviction is dropped | 3 tests |
| eviction at the region's end instead of its first base | 8 tests |
| build over the analysed region instead of the building region | 8 tests |
| the failed spans are not gathered | 2 tests |
| the held count sums only the first sample | `the_held_count_sums_every_samples_window` |
| a one-base gap between consecutive building regions | 8 tests |
| the division does not order an inverted region's ends | `the_division_reads_an_inverted_region_as_the_ground_it_names` |
| `min_alt_obs` replaced by its default | 2 tests |

**One thing is deliberately not pinned**, and the test says so: evicting *after* the cover
instead of before gives the same output and the same final window. What it costs is the sweep —
measured, the pre-cover eviction takes the window from three records to one at every region.

**One fixture defect, found by running it:** the first version gave the dotted sample reference
base `G` where the deletion's reference was `A`, and two of its positions sit inside the
deletion's span. `build_region` refuses that — two members of one locus disagreeing on the
reference means the samples were called against different references — so the fixture failed in
*both* drivers. The whole fixture's reference is now `A`.

## 9. Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` → `164 passed; 0 failed`. `--all-targets` clippy is red
on this branch for 49 pre-existing reasons in other modules — the standing item.
