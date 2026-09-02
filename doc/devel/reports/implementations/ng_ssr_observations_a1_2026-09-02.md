# ng STR observations — A1: the merge carries the locus kind

*2026-09-02. Step A1 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §4](../../ng/spec/run_ssr_observations.md). Branch `ng-ssr-observations`.*

## Plan

The merge closes repeat-tract loci deliberately — the width bound exempts them, and the walk
asserts that no locus mixes a tract with ordinary sequence — and then drops what kind of
ground it was at the moment it assembles the cohort's evidence. With the kind went the
motif, so the period, so the repeat counts, and a caller handed a cohort observation could
not tell a tract from a stretch of ordinary sequence.

Two fields, one clone:

- `ClosedLocus<'a>` gains `pub kind: &'a LocusKind`, borrowed from the observation that
  opened the locus — the same value `judge` was already handed.
- `CohortObservation` gains `pub kind: LocusKind`, cloned from the closed locus in
  `CohortObservation::over`.

## Assumptions

**The spec says "`ClosedLocus` exposes the locus kind" without saying owned or borrowed;
this borrows.** The plan's own wording settles it — *"the clone is per built locus"* — and
the arithmetic agrees: a closed locus is produced for every locus the walk closes, the
refused ones included, while `CohortObservation::over` runs only on the loci the caller
undertakes to build. Borrowing on `ClosedLocus` costs nothing at a refused locus and the
type already names a lifetime, so no signature gained one.

**The hand-built fixture helper in `build.rs` takes the kind from the first observation its
members hold, which is not the walk's rule.** The walk takes it from whichever observation
starts earliest and can afford to, because it asserts the members agree by discriminant. A
fixture can present members that disagree; this helper then answers with the first sample's
rather than refusing. Nothing tests that case and nothing needs to — `close.rs`'s
`a_locus_mixing_an_str_tract_and_a_generic_observation_is_refused` is what says it cannot
arise from real input. Written into the helper's comment rather than left to be rediscovered.

## Changes made

| file | change |
|---|---|
| `src/ng/run/cohort_merge/close.rs` | `ClosedLocus::kind: &'a LocusKind`, set from the opening observation |
| `src/ng/run/cohort_merge/build.rs` | `CohortObservation::kind: LocusKind`, cloned in `over`; fixture helper carries a kind |
| `src/ng/calling/allele_candidates/mod.rs`, `src/ng/calling/evidence_shaping.rs`, `src/ng/run/cohort_merge/{organise,serial}.rs`, `src/ng/run/records/tests.rs` | test fixtures name `LocusKind::Generic` |
| `tests/ng_calling_loop_calls_genotypes.rs`, `tests/ng_candidate_selection_truth_recall.rs` | the same, at the integration fixtures |
| `examples/ng_candidate_selection_probe.rs` | the probe's sample-restricting copy carries the kind through |

**One file outside this plan's ownership table was touched**:
`examples/ng_candidate_selection_probe.rs`, which the parallel loop plan owns. The edit is
one line forced by the new field — a struct literal that no longer compiles without it — and
carries the kind rather than defaulting it.

No verdict moves and no behaviour changes: nothing yet reads either field.

## Tests added

- `close.rs::a_closed_locus_carries_the_kind_its_members_share` — a tract locus and a generic
  locus through the walk; the tract arm compares the **whole** `LocusKind`, so a kind rebuilt
  as a bare `Ssr` with an empty motif fails it where a discriminant comparison would pass.
- `close.rs::a_refused_locus_still_states_its_kind` — a tract nobody varied at comes out
  `TooQuiet` and still carries its motif, which is what makes the borrow's argument checkable:
  the field is populated on loci `CohortObservation::over` never sees.
- `build.rs::the_assembled_locus_carries_the_kind_the_generator_minted` — the motif and both
  flanks survive `CohortObservation::over`, against a generic control assembled from the same
  region.

## Validation

All run in the dev container (`./scripts/dev.sh`).

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --tests --examples --all-features --no-fail-fast` — 5,915 passed, 0
  failed, 14 ignored in the library suite; every integration target green; **two example
  targets red before this step and after it** (below).
- The three new tests by name: 3 passed.

**Two pre-existing failures, neither this step's and neither fixed here.**

`examples/ng_generic_loci_dump.rs` (11 of 13 tests) and `examples/ng_ssr_loci_dump.rs` (all
10) fail at fixture setup with
`RepeatCatalogError::NotFound { path: ".../ref.fa.repeats.parquet" }`. Their `fixture` helper
writes a FASTA and a BAM into a temporary directory and never builds a repeat catalog beside
them, which stopped being optional at `8ffb0f84` (*"the live tandem-repeat scan is gone"*,
2026-08-11): nothing scans a reference for tandem repeats at call time any more, so a dump
over a reference with no catalog file has nowhere to read them from. Bisected — `bb27dd42`
green, `8ffb0f84` red, and both targets red at `main`'s tip (`e270ce14`) without this step's
commit.

**One of the two is spec §10's own oracle for this milestone**: *"a tract cohort observation
reaching a probe carries the same motif the generator minted, on real data (`ng_ssr_loci_dump`
ground truth beside the merge's output)"*. That check cannot be run until the fixture builds a
catalog. The unit tests added here pin the field on synthetic loci; the real-data
confirmation is owed.
**And a third, in the benches.** `cargo test --all-targets` also builds the criterion benches,
and `benches/psp_writer_perf.rs`
panics at line 386 — `index out of bounds: the len is 3300000 but the index is 3300000` — in
`psp_writer_phases/flush_block_one`, which walks its fixture until the block fills and runs
off the end when it does not. Confirmed pre-existing by stashing this step's edits and
re-running `cargo test --bench psp_writer_perf` on the untouched tree: the same panic, at the
same line. The bench is over production's `.psp` writer and names nothing this step touched.

`cargo doc --no-deps` reports 26 unresolved intra-doc links, all pre-existing and none in
`close.rs` or `build.rs`; the links added by this step resolve.

## Tradeoffs and follow-ups

- **Nothing reads the field yet.** The driver still sends every cohort observation down the
  SNP/indel path; the branch on the kind is the loop plan's Milestone C, and this plan's C2
  puts a counted set-aside in front of it in the meantime.
- **`LocusKind::SsrBundle` carries no payload**, so a bundle observation states its kind and
  nothing else. Deliberate — spec §8 defers the bundle payload to its own document.
