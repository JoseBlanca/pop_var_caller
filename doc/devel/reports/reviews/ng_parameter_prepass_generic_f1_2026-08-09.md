# F1 review — synthesis and fixes applied

**Date:** 2026-08-09. **Reviewed:** `257ccdaa`. **Two agents in isolated worktrees**: mutation
testing, and design/contract. Verdicts: Approve-with-changes ×2. **2 Blockers, 3 Majors**,
plus one wrong claim of the author's.

---

## 1. The two Blockers, and both are the fixture choosing what it could see

**The whole `Fitted` inbreeding arm was dead to the suite.** All four F1 tests built their
config with `InbreedingMode::Supplied`, to keep the runs model out of the way. Replacing the
entire `InbreedingMode::Fitted => …` arm with `unreachable!()` left **297 tests green**, as did
suppressing the runs model with a bare `Ok((None, None))`. A build that silently reported no
`F` for every fitted sample would have shipped.

This is the same fault as E2's twenty-reads world: the fixture was chosen to make the test easy
to write, and it could not see the thing the step is. Closed by
`a_fitted_coefficient_is_attempted_rather_than_reported_absent` — the fixture's 10,001 one-base
sites lie in a single 100 kb window against a floor of 3,000, so the runs model refuses, and
*reaching that refusal* proves the arm ran, the ploidy-2 lookup found its entry and
`library_noise` was built.

**`covered_positions` could be pointed at the wrong table with nothing to notice.** Summing the
read-group table instead of the windowed one keeps the suite green — which is exactly the error
the method's own doc says the two tables exist to prevent. Every fixture gave each site one
read group, including the three-group test whose groups cover *disjoint* thirds, so the two
tables agreed everywhere the suite looked. Closed by
`a_site_two_libraries_covered_counts_one_position_and_not_two`, the only new fixture shape the
review needed: 20,002 against the correct 10,001.

## 2. Three Majors

**The runs model was handed library shares pooled across ploidies.** `library_shares` sums
`total_reads()` over `(group, ploidy)` ignoring the ploidy, while the runs model walks the
diploid windows and nothing else. Measured on a fixture whose two arms come from different
libraries: 0.5/0.5 shares where the diploid arm's reads are entirely one of them, which at
Phred 20 against Phred 30 puts the share-weighted rate **5.5× above the truth** — on the one
model whose job is separating real heterozygotes from error. Fixed by `library_shares_over`,
which takes a predicate on ploidy.

**`InbreedingMode::Supplied` overflowed the windowed table's `u32` counters on any real
genome.** Reproduced: a panic at **143,165,576 sites at depth 30 in one cell**. Pre-existing
from Milestone C, and F1 is what made the mode reachable from production. Fixed on the owner's
call with a separate `u64` table keyed by ploidy alone — see §4 for what it nearly broke.

**Nothing read `coupled_fit`'s termination.** Reporting a constant
`FitTermination { iterations: 999, converged: false }` left the suite green: the equality test
between the two entry points cannot see it, because both sides carry the same constant. A fit
that ran out of iterations would be reported as converged. Closed, along with
`accumulators()`'s edges-sharing contract, which had nothing executable behind it — no test
merged two shards, which is the whole reason the second entry point exists.

## 3. A wrong claim of the author's

`LocusGeneration` carried a rendered `String` rather than the walk's error, on the stated
grounds that `ParameterEstimationError`'s `Clone` and `PartialEq` were relied on by its tests.
**They are not.** A reviewer measured it and so did the author: removing both derives compiles
with **zero** errors across `--all-targets`. The variant now carries
`#[source] LocusGenerationError`, which matters because five of its six variants name the
`GenomeRegion` where the walk broke — reachable before only by parsing prose.

## 4. What the cheap fix nearly broke, and why the review's suggestion was not taken verbatim

The owner chose the separate-`u64`-table fix over widening `by_window`, because the collapsed
arm holds one table per ploidy — a few kB — against `by_window`'s ~8,000, which is the 37 MB
`Supplied` mode exists to save.

**But `fit_coupled` derived the set of ploidies from `windowed_histograms().keys()`.** With the
collapsed sites moved to their own map, that is empty in `Supplied` mode — so the coupled fit
would have found no ploidies and fitted nothing at all, *silently*, since an empty map is not
an error. `GenericAccumulators::ploidies()` answers from whichever map the mode uses. Without
it the fix would have traded a loud panic for a silent no-op, which is worse.

The existing collapse test asserted the old structure. Rather than relocate its assertion it
now asserts the property the fix is *for*: the `u32` map must receive **nothing** in this mode,
because a site landing there is a site on the path that overflowed.

## 5. One reviewer suggestion recorded rather than applied

**The `library_noise` transposition survives.** Pairing shares with rates in reverse is a
genuine transposition at two or more libraries, and it survives even with a two-library fitted
test, because `fit_inbreeding` refuses on the window count before it touches the noise.
Consuming that noise needs a fixture spanning 3,000 windows — 300 Mb — so it belongs to F2 or
to a harness, not to F1. Recorded so it is not assumed covered: this is consistent with the
module's own claim that no identity of the scoring rule can catch a transposition, and it means
the keyed construction the doc calls "the only thing that prevents it" has nothing checking it.

## 6. Two things for the owner

- **`config.read_groups` made no observable difference on the reviewer's probe.** Two read
  groups on the same 10,001 sites, one supplying every alternative read and the other none —
  the most asymmetric two-library shape available — returned `GenericSampleParameters` that
  compare **equal** whether `multi_library` is true or false. Either the attributed arm changes
  nothing measurable at that fixture size, or the fixture does not reach it. Worth settling
  before F2 chooses its fixtures.
- **An empty sample panics rather than erroring.** `fit_coupled` asserts the sample has a read
  group with reads. Documented and intentional, but in a cohort run a sample whose requested
  regions yield no generic loci aborts the process rather than failing that one sample. Worth
  an error variant when the cohort gather lands.

## 7. Design-doc drift, added to the owner's list

- **arch §1.1** mandates that a `LocusGenerationError` "is fatal and propagates" without saying
  as what; `ParameterEstimationError::LocusGeneration` now exists and is not in the doc.
- **arch §5.4's error enum block** omits four live variants — `GenotypeFrequenciesOffSimplex`,
  `InbreedingStartsDisagree`, `LocusGeneration` — shows two others without the `floor` field
  they carry, and shows a `Domain(#[from] DomainError)` that is now a struct variant.
- **arch §3's accumulator surface** owes rows for `covered_positions`, `ploidies`,
  `read_group_histograms`, `windowed_histograms` and `inbreeding_mode`, and a correction to
  `adjustments`' return type.
- **arch's module table** owes a `generic/estimate.rs` row.
- **arch §5.4's "supply it or accept the failure"** now has a stated consequence: the failure is
  the whole sample's, decided by the owner on 2026-08-09 because `F` is a prior the calling step
  needs.

## 8. Verification

Every fix was re-run. F1's tests: 4 → **8**. Suite 3,199 → 3,203. Gates: `fmt` clean, `clippy
--all-targets --all-features -D warnings` clean, `doc --no-deps --lib` at its 12-unresolved-link
baseline with none in this module.
