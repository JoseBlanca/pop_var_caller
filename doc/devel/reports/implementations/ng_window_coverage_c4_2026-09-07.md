# window coverage — C4: the pair on the locus, and beside every record the run writes

**Date:** 2026-09-07
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone C, step C4
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.5, §10
**Review:** [ng_window_coverage_c4_2026-09-07.md](../reviews/ng_window_coverage_c4_2026-09-07.md)
**Branch:** `ng-window-coverage`
**Builds on:** [C2](ng_window_coverage_c2_2026-09-06.md), [C3](ng_window_coverage_c3_2026-09-06.md)

## The answer

**Every record the run writes now carries, per sample, the mean read depth and GC of the 500-base
window centred on the locus's first base — and it carries them *beside* the record, not in it, so
the VCF is unchanged.** Measured on the six tomato accessions over two 100 kb intervals:

- **The two modes agree bit for bit**: 13,866 rows — 2,311 records × 6 samples — identical between
  a run over the alignment files and a run over the stored ones. This is what the plan asks C4 to
  prove, and the mode-equivalence oracle now makes the comparison itself.
- **Every one of those 13,866 sample-loci either matches the whole-store recomputation bit for
  bit, or has nothing to report on both sides**: 13,589 carry a window and agree, 277 have none on
  either side, **0 disagree, and 0 are ones the run has no window for where the walk does**.
- **The VCF is unmoved**: 2,311 records, sha256
  `84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d` on both routes, as at every
  step since A1.

**The twelve the run had no window for at C3 are gone, and they were never written records.** C3's
comparison was over the 26,754 sample-loci the merge *builds*; this one is over the 13,866 it
*writes*. The two positions those twelve sat at — the last two loci of the last analysed interval
— establish no variant and become no record. So at the seam the filter will actually read, the
run's windows are complete on this data.

## Where the pair lives, and why in three shapes

The measurement crosses three boundaries and takes a different shape at each, because each side
indexes differently:

| where | shape | indexed by |
|---|---|---|
| the merge's cohort locus | `Vec<WindowCoverage>` on `CohortObservation` | the **covering** samples, parallel to `per_sample` |
| the evidence gathered for output | one `WindowCoverage` on `SampleEvidenceForOutput` | the **run's** samples, dense |
| beside the record at the sink | `&[WindowCoverage]` | the run's samples, dense |

**The sparse-to-dense translation is the one place the two orders meet**, and it happens where the
same loop already translates the read counts (`records.rs`, `evidence_for_output`): the merge's
entry *i* names run sample `support.sample`, and the window at *i* goes to that row. A sample that
covered nothing at the locus has no entry in the first and keeps an absent pair in the second,
which says the same thing.

**The record itself gains nothing**, which is what keeps the VCF byte-identical. The pair travels
beside it, through the sink's second parameter, so the hidden-duplication filter has what its
spill needs without re-deriving anything from the record; the two CLI adaptors ignore the slice,
and say so.

## The decision this step had to make, and how it was settled

**`WindowedCohort::window_coverage_at` answers `None` for two different things** — "this sample has
no window here", which is ordinary, and "this view carries no measurement at all", which is a gap.
Three construction sites pass no measurement, and one of them is `merge_cohort_serially`, the
in-memory driver every cached-driver test is compared against. Left alone, the drivers' agreement
comparison would have **passed by comparing absence against absence**.

**Settled by refusing to compare it rather than by splitting the type.** Spec §3.5's type block is
`Vec<WindowCoverage>` and the consumer is being built on another branch, so widening it to an
`Option` would ripple past this step for a distinction only the tests need. Instead:

- `cohort_merge::fixtures::render` — the one place that decides what "the same answer" means
  across this module's drivers — now **destructures `CohortObservation`** and leaves the window
  coverage out, so a field added later is a compile error here rather than a comparison that
  quietly stops covering it;
- and `serial::tests::the_cached_driver_measures_windows_where_the_in_memory_one_has_none` says
  the asymmetry out loud: the cache's pairs are numbers, the in-memory driver's are all absent,
  both are one entry per covering sample, and **some locus carries a number** — which is what stops
  it passing over an outcome of absences.

## Assumptions and deviations, all minor and all recorded

1. **`CohortObservation::over` takes the window it was resolved against**, rather than the field
   being filled after construction. Two-phase construction is a field a later caller can forget;
   the cost is a `&no_measurement()` argument at about thirty fixtures, which is one shared
   helper in the module's test scaffolding.
2. **`SampleEvidenceForOutput` drops `Eq` from its derive**, because a derive cannot reach past a
   field whose type withholds it — and `WindowCoverage` withholds it deliberately, since bitwise
   equality separates `+0.0` from `-0.0`. `PartialEq` still means what it meant: `WindowCoverage`'s
   own compares by bit pattern, which is what makes two absent windows equal. Nothing in the tree
   compares two of these, and the only type that holds them never derived `Eq` either.
3. **The recorder moved from the locus builder to the record sink**, which C3's design review had
   recommended for exactly this step: spec §10 says "at every written record", and this is where
   that is decided. It takes the dense slice now rather than a `WindowedCohort`.
4. **The whole-store comparison collapses "absent" and "missing" into one answer.** By the time a
   window reaches a record, the merge's three-way distinction is one pair of `NaN`s
   (`CohortObservation::window_coverage`'s own doc says why), while the recomputation still answers
   `None` where it finalised no centre. Comparing the spellings rather than the fact reported a
   disagreement at **277** sample-loci on the first run — every one of them `run=absent` against
   `walk=none` — which is a difference in how nothing is written down, not in what was measured.

## Changes made

- **[`build.rs`](../../../../src/ng/run/cohort_merge/build.rs)** — the field, and `over` reading it
  at the locus's first base.
- **[`cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)** — `render` destructures
  and excludes.
- **[`serial.rs`](../../../../src/ng/run/cohort_merge/serial.rs)** — the test that says what
  `render` leaves out.
- **[`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)** — the
  builders' differential compares through `render` rather than `Debug` on the whole outcome.
- **[`recorded_windows.rs`](../../../../src/ng/run/cohort_merge/recorded_windows.rs)** — the
  recorder over the dense slice, called from the record sink rather than the locus builder.
- **[`assemble.rs`](../../../../src/ng/vcf/assemble.rs)**,
  **[`records.rs`](../../../../src/ng/run/records.rs)** — the dense pair and its translation.
- **[`callers.rs`](../../../../src/ng/run/callers.rs)**,
  **[`psp_caller.rs`](../../../../src/ng/run/psp_caller.rs)**, both CLI commands — the sink's
  signature and the scratch buffer that fills it.
- **[`ng_mode_equivalence_oracle.sh`](../../../../scripts/ng_mode_equivalence_oracle.sh)** — both
  routes record their windows and the two files are compared row for row, refusing an empty
  comparison **and one in which no row carries a number**.
- **[`ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)** — the
  comparison collapses the spellings of "nothing to report".

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched.

## What the review changed

Two agents; the full report is [beside this one](../reviews/ng_window_coverage_c4_2026-09-07.md).
**One Blocker, found by both, four Major, twelve Minor** — all applied. The three that were
defects rather than clarity:

- **Two integration tests stopped compiling**, and `cargo test --lib` does not build them, so
  twenty assertions had quietly stopped running. `cargo check --lib --tests` is part of what green
  means here now, and `WindowedCohort::with_nothing_measured()` gives the out-of-crate callers one
  documented way to say "no measurement" instead of a struct literal per file.
- **The extended mode-equivalence oracle passed with the measurement entirely off.** Two files of
  absent pairs are byte-identical, so it agreed loudest exactly where the measurement was missing —
  the same "absence against absence" failure this step's own design decision was taken to avoid,
  one level up. It now requires a row that carries a number, and says how many: `13,866 rows,
  13,589 of them a measurement`, which cross-checks against the whole-store comparison's own count.
- **"The locus's first base" had no test.** Every fixture with a real measurement uses a one-base
  locus, where first and last coincide, so reading the *last* base passed the whole library suite —
  and the mode-equivalence oracle cannot see it either, because both modes would read the same
  wrong position. A five-base locus with a different window at each end now pins it.

And two smaller ones worth carrying: the merge's two per-sample lists meet at exactly one index and
only one side was checked, so a locus whose lists came apart failed as a bare bounds panic blaming
the reader; and three doc comments had stopped describing the code — a `-` row the writer can no
longer produce, an `Eq` removal blamed on `NaN` when `NaN` is the one thing the comparison handles,
and `render`'s own doc arguing against the code beneath it.

**The mutation table is the useful part of the report.** Three of the four defects that survive the
library suite are caught only by the real-data probe, and none by the mode-equivalence oracle —
which compares two routes that share the defect. The oracle proves the two modes agree; the probe
proves what they agree on is right.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,373 passed, 0 failed, 15 ignored** (6,370 before this
  step); `cargo test --all-features --example ng_window_coverage_probe` — 7 passed.
- `cargo check --lib --tests --all-features` — clean. **New to this step's checklist**, because
  the review found it failing and `--lib` cannot see it.
- `cargo test --test ng_calling_loop_calls_genotypes --test ng_candidate_selection_truth_recall` —
  19 passed, 1 failed, and the failure is the pre-existing
  `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`.
- `cargo clippy --lib --all-features --bins --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check` — every file this step touched is at the hunk count it had before it.
- **The standing oracle**: 2,311 records, sha256 `84ad19c2…0590d` on both routes, **and the two
  modes' window coverage identical bit for bit over 13,866 rows, 13,589 of them a measurement**.
- **The whole-store recomputation**: 13,589 of 13,866 sample-loci agree bit for bit, 277 have
  nothing to report on either side, 0 disagree.

## Tradeoffs and follow-ups

- **Nothing reads the pair yet.** The hidden-duplication filter is its consumer and is built on its
  own branch; until then both CLI adaptors take the slice and drop it.
- **The histograms are still inside the cache.** C5 takes them out, and is also what closes the
  last half-window of each sample's stream — the centres no cover can finalise.
- **The recorder writes one row per sample per written record**, which at 63 accessions is 63 rows
  a record. It is opt-in and no run that is not measuring pays for it.
