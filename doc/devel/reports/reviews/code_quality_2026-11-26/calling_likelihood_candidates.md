# Code-quality map: read likelihood, candidate selection, genotype prior, site quality

Area: `src/calling/likelihood/`, `src/calling/allele_candidates/`,
`src/calling/genotype_prior/`, `src/calling/quality/` — about 38,000 lines, of which roughly
15,000 are production code and the rest tests. Read-only review, 2026-11-26, from the code
alone (nothing was built or run).

## Summary

1. **The evidence hand-off types live in the run driver, so the library imports `crate::run`.**
   `CohortObservation` and the five types it is made of are defined in
   `run/cohort_merge/build.rs`, and the merge's support thresholds (`MinAltReads` and its two
   parts) in `run/cohort_merge/mod.rs`. Nine non-test imports across this area, plus two in
   `calling/` outside it, point from the calling library into `run`, while `run` drives
   `calling` — a two-way dependency. Moving the plain data types to a new top-level module
   (`src/cohort_observation.rs`) removes every one of them and changes no output. This is
   the one high-value finding.
2. **The cap-and-admit tail of candidate selection is written twice** — about 55 lines,
   once in `generic.rs` and once in `ssr.rs`, differing only in where the leftover rows and
   the locus kind come from. One shared helper would hold the rule that decides the `ALT`
   column order in one place.
3. **The partial-read fit test is written twice**: the STR selection's `partial_read_fits`
   is the same comparison as the generic likelihood's `carrier_sequence_starts_with` /
   `carrier_sequence_ends_with`.
4. **No confirmed bug.** I traced selection, the site-quality fold, the binomial tails and
   the contamination leave-one-out end to end and found the edge cases guarded. Two test
   notes are under Bugs.
5. **Readability: stale references** — a deleted `hardy_weinberg` module described as
   present, three `examples/*.rs` programs that no longer exist, two nonexistent methods
   linked from a doc comment, a test cited by the wrong name, and production-era file paths
   (`src/var_calling/…`, `src/ssr/…`) that were deleted on 2026-09-13.
6. **Small clean-ups**: two public functions that only tests call (one of them re-exported
   at `calling::`), a numeric helper (`log_sum_exp_2`) living inside the prior's model file,
   and a flank-respelling helper that selection borrows from the likelihood's internals.

---

## Bugs

No finding here meets the brief's bar (a concrete input that gives a wrong call, a crash or
lost data). What I checked, so the next reader need not repeat it:

- `allele_candidates`: the fold, cap, admission and leftover (`mod.rs:893-1194`,
  `generic.rs:81-192`, `ssr.rs:1087-1291`). The scratch is reset per locus by
  `summarise_alleles` → `reset_for`, so `ranked_table_indices.extend` in `generic.rs:124`
  starts from empty. Ordering is total (bases break the last tie), so output does not depend
  on sample order or thread count.
- `quality/mod.rs:451-863`: the site-quality fold. The band bookkeeping and stale-entry
  clearing are correct; `peak` is never zero because each sample's taps are normalised to a
  maximum of one; an all-`−∞` sample is caught by an assertion with a named cause.
- `quality/artifact_correction.rs:144-300`: the two-sided binomial tail. Both binary
  searches keep their invariant, and the end clamps keep the Phred finite.
- `likelihood/mod.rs:717-1052`: the contamination tables. Leave-one-out subtraction is
  floored at zero; a non-finite batch total is refused.

Two test notes (in scope per the brief: tests that may not test what they claim):

### T1. The generic row's tests mostly go through a wrapper production never calls

- **Where** — `src/calling/likelihood/generic.rs:771` (`genotype_log_likelihood_row`); the
  production path is `fill_generic_emissions` + `assemble_genotype_log_likelihood_row`, called
  from `src/calling/inference/summarise_condition.rs:1867` and `:2354`.
- **What** — the wrapper passes `ReferenceBesideLocus::NONE` (no reference beside the locus),
  and production passes the real flanks. For a partial read whose bases run past an allele,
  `carrier_sequence_starts_with` with an empty flank says *incompatible*, and with the real
  flank it can say *compatible*. So any test that scores a partial read through the wrapper
  (seven test call sites: `generic.rs:1803, 2391, 2716, 2770, 2896, 3710, 4493`; I did not
  check which of them carry partials) exercises a stricter comparison than the caller runs. Dedicated tests do cover the flank case
  (`after_the_giab_deletion`, `generic.rs:4165`), so this is a gap in emphasis, not missing
  coverage.
- **Proposal** — mark the wrapper `#[cfg(test)]`, drop its re-export at
  `src/calling/mod.rs:79`, and state in its doc comment that it scores without flanks.
- **Downside** — none of substance.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### T2. A "parity with production" test now pins a literal against a literal

- **Where** — `src/calling/likelihood/mod.rs:2932` (`the_floor_still_equals_productions`).
- **What** — it compares `MIN_BASE_ERROR` against a local `const` set to the same `1e-12`.
  The production code it once tracked was deleted on 2026-09-13, so the test can only fail if
  someone edits the constant, which is a deliberate act and not drift.
- **Proposal** — leave it or delete it; if kept, rename it to say what it does (pins the
  value). Its doc reference at `mod.rs:204` also has the name wrong (`the_floors_still_equal_productions`).
- **Effect** — none. **Cost** — S. **Confidence** — high.

---

## Architecture

### A1. Move the merge-to-calling hand-off types out of `crate::run` (the known lead)

- **Where**
  - Types: `src/run/cohort_merge/build.rs:1289` `CohortObservation`, `:1434` `SampleSupport`,
    `:1559` `SupportedAllele`, `:1605` `PartialObservation`, `:1683` `LocusBorder`, `:1751`
    `AlleleSupport`; `src/run/cohort_merge/mod.rs:574` `MinAltObs`, `:627` `MinAltReadShare`,
    `:725` `MinAltReads` (plus `DEFAULT_MIN_ALT_OBS`, `DEFAULT_MIN_ALT_READ_SHARE`).
  - Imports in this area: `allele_candidates/generic.rs:15`, `allele_candidates/mod.rs:44-46`,
    `allele_candidates/ssr.rs:46-47`, `allele_candidates/deciding_partial_reads.rs:47`,
    `likelihood/generic.rs:13`, `likelihood/mod.rs:142`. Outside it, in `calling/`:
    `evidence_shaping.rs:77`, `inference/mod.rs:67`.
- **What** — these are the *contract* between the cohort merge (producer) and calling
  (consumer): plain data with public fields and no behaviour the calling side needs except
  `flank_deletion_of`. Because they are defined inside the run driver, the calling library
  depends on `run`, and `run` depends on `calling` (`run/callers.rs:33` imports the selectors).
  Nothing in the type definitions needs `run`: their fields use only `crate::types`
  (`GenomeRegion`, `ReadGroupId`), `crate::locus_generation` (`LocusKind`, `FlankDeletion`,
  `WitnessedLocusPositions`, `HeldFlank`) and `crate::window_coverage`.
- **Where it should live — the options weighed**
  1. **A new top-level module, `src/cohort_observation.rs` (recommended).** Holds the six
     evidence types and the three threshold types. Producer (`run::cohort_merge`) and
     consumers (`calling`, `vcf`, `psp` docs) both depend on it; it depends only on
     `types`, `locus_generation` and `window_coverage`. Its name says what it holds.
  2. `crate::locus_generation`. It already holds the per-sample evidence types
     (`SequenceObservation`, `LocusKind`, …), so the move is natural in kind, but that module
     is about *one sample's* reads; putting the cohort-level type there blurs it, and the
     module is already large.
  3. Inside `crate::calling` (consumer owns the contract). Correct dependency direction
     (`run` → `calling` already exists), but the merge's own type then lives in the calling
     library, and the merge's tests in `run` would import from `calling` to build fixtures.
  4. Move all of `cohort_merge` out of `run`. Not now: it imports `run::explain`,
     `run::side_pool`, and `cli::calling_run::round_shape_for`, so it is coupled to the run
     driver in ways the type move does not need to touch.
- **How** — move the struct and enum *definitions* only. The constructors that need run-side
  types (`CohortObservation::over`, `over_from_cursors`, `over_looking_up`, which take
  `ClosedLocus` and `WindowedCohort`) stay in `build.rs` as an inherent `impl
  CohortObservation` block — Rust allows inherent impls anywhere in the defining crate — so the
  producer code does not change. `build.rs` re-exports the moved names for one step if that
  keeps the diff small, then the re-export is removed.
- **Downside / why it might be wrong** — about 520 lines move, most of them doc comments, and
  those docs link to `build.rs`-private items (`partials_of_sample`, `alleles_of_sample`,
  `MemberPlacement`, `SampleMembers`, `AlleleTable`) that would become plain-text references.
  The docs also mix contract (what a field means) with producer detail (how the merge fills
  it); the move is a chance to split them but makes the change larger if done at once. The
  threshold types (`MinAltReads`) are less obviously "evidence"; they could stay in
  `cohort_merge` with only the evidence types moving, but then `allele_candidates` and
  `inference` keep a `run` import, so I recommend moving them together.
- **Effect** — none on output, speed or memory (definitions only). **Cost** — M.
  **Confidence** — high.

### A2. Selection borrows a flank-respelling helper from the likelihood's internals

- **Where** — defined `src/calling/likelihood/ssr_emission.rs:115` (`respelled_over`), used by
  `src/calling/allele_candidates/ssr.rs:43` (in `one_base_apart`, `:1399`) and
  `src/calling/likelihood/ssr.rs`.
- **What** — `respelled_over` rewrites a sequence over a wider flank-deletion extent. It is an
  operation on `FlankDeletion`, not on the emission model, but selection reaches into
  `likelihood::ssr_emission` for it, which makes candidate selection depend on the STR
  emission module.
- **Proposal** — move it next to `FlankDeletion` in `src/locus_generation/mod.rs` (as a method
  or a free function).
- **Downside** — one more item in an already large `locus_generation/mod.rs`; the coupling it
  removes is one function.
- **Effect** — none. **Cost** — S. **Confidence** — medium.

### A3. `log_sum_exp_2` lives inside the prior's model file and is imported by site quality

- **Where** — `src/calling/genotype_prior/dirichlet_multinomial.rs:361`; imported at
  `src/calling/quality/mod.rs:37`. Site quality also has its own `log_sum_exp_over`
  (`quality/mod.rs:856`).
- **What** — two general numeric helpers living in domain files; site quality depends on the
  genotype-prior module only to reach one of them.
- **Proposal** — move both, bodies unchanged, into `src/float.rs`, which already holds the
  crate's `ln`/`exp`. **Do not** merge them with `paralog/locus_score.rs:700`
  (`log_add_exp`): that one swaps `log1p` for a cubic series and would change output.
- **Downside** — `float.rs` is today a thin wrapper over libm; this widens its role slightly.
- **Effect** — none if bodies are moved verbatim. **Cost** — S. **Confidence** — high.

### A4. Looked at: the dependency between `likelihood` and `genotype_prior`

`likelihood/mod.rs:134` imports `COUNT_PATH_DESYNC_THRESHOLD` and `SampleAlleleCopies` from
`genotype_prior`. One direction only, and both items are about the same expected-copy counts
the prior owns. Leave it.

---

## Duplication

### D1. The cap, the base-order sort, admission and the leftover are written twice

- **Where** — `src/calling/allele_candidates/generic.rs:136-191` and
  `src/calling/allele_candidates/ssr.rs:1229-1290`.
- **What** — about 55 lines each, the same sequence: compare survivors to the cap, rank with
  `compare_best_first`, record the cut indices, truncate, sort survivors by bases, admit them
  into `CandidateAlleles` with a remap, then compute each sample's leftover. The differences
  are three parameters: the locus kind (`LocusKind::Generic` vs `ground.locus_kind()`), the
  rows the leftover reads the read count from (the same rows vs `merge_rows`), and `admit` vs
  `admit_with_flank_deletion` — which are identical when the flank deletion is `NONE`, as it
  always is on the generic path (`src/calling/mod.rs:190-206`, `:279-288`). This tail decides
  the `ALT` column order and which samples are no-called, so two copies are two places a fix
  has to land.
- **Proposal** — one `pub(super) fn cap_and_admit(observation, merge_rows, locus_kind,
  config: &CandidateSelectionConfig, scratch) -> LocusSelection` in `allele_candidates/mod.rs`,
  called by both paths after they have filled `ranked_table_indices`.
- **Downside** — the module doc says the two paths are deliberately two functions rather than
  one trait. This keeps that: each path still chooses *which* sequences are candidates; only
  the shared tail moves. A future STR-only rule in the tail would need a parameter.
- **Effect** — none (same operations, same order). **Cost** — S–M. **Confidence** — high.

### D2. The partial-read fit test exists twice

- **Where** — `src/calling/allele_candidates/deciding_partial_reads.rs:207`
  (`partial_read_fits`) and `src/calling/likelihood/generic.rs:608`, `:617`
  (`carrier_sequence_starts_with`, `carrier_sequence_ends_with`).
- **What** — `partial_read_fits` with `HeldFlank::Left` is exactly
  `carrier_sequence_starts_with(allele, far_flank, read)`, and with `HeldFlank::Right`
  exactly `carrier_sequence_ends_with(far_flank, allele, read)`: both take the shorter of read
  and allele, require agreement over it, and require the rest of the read to match the flank.
  Traced both branches, including a read shorter than the allele.
- **Proposal** — have `partial_read_fits` compute `far_flank` as it does now and call the two
  generic functions (they are already `pub(crate)`).
- **Downside** — couples STR selection to the generic likelihood file; if A2 is done, move the
  two helpers to `locus_generation` alongside `respelled_over` and both callers import from
  there.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### D3. Leave-one-out subtraction written three times

- **Where** — `src/calling/genotype_prior/mod.rs:734` (`fill_sample_concentration`),
  `src/calling/likelihood/mod.rs:832` (`fill_contaminant_allele_frequencies`), `:988`
  (`fill_own_batch_contaminant_row`).
- **What** — each subtracts a sample's own expected copies from a total, checks the result
  against `COUNT_PATH_DESYNC_THRESHOLD` in a `debug_assert!`, and floors it at zero. Three
  lines of body plus a seven-line assertion message each.
- **Proposal** — optional: one `leave_one_out(total, own, allele) -> f64` beside
  `COUNT_PATH_DESYNC_THRESHOLD`. Better first step: see R2, which removes one of the three.
- **Downside** — the bodies are short; a helper saves little and ties the prior to the
  contamination code. **Leave it** unless R2 is done and two copies still read as one idea.
- **Effect** — none. **Cost** — S. **Confidence** — medium.

### D4. "Sum one allele's reads over its read-group rows" written four times

- **Where** — `src/calling/allele_candidates/mod.rs:915`, `:941`, `:1089`; and
  `src/calling/allele_candidates/ssr.rs:1382` (`a_minority_beside_it_everywhere`).
- **What** — the same saturating fold over `row.support.num_reads`; the `ssr.rs` copy uses a
  non-saturating `.sum()` over a filter instead.
- **Proposal** — one `fn pooled_reads(rows: &[SupportedAllele]) -> u32` in `mod.rs`.
- **Downside** — trivial gain. Do it only alongside D1.
- **Effect** — none (overflow needs four billion reads). **Cost** — S. **Confidence** — high.

---

## Readability

### R1. Stale references to deleted code

- **Where and what**
  - `src/calling/genotype_prior/mod.rs:72-76` lists a `hardy_weinberg` module ("the
    comparator … `PlugInWrightPrior`"); `:109-112` counts it among four descendant files;
    `dirichlet_multinomial.rs:352-354` says `log_sum_exp_2` has three callers including it.
    There is no such file in `src/calling/genotype_prior/`.
  - `src/calling/allele_candidates/mod.rs:739-744` links
    `SelectionScratch::best_within_sample_share_of` and `SelectionScratch::cohort_reads_of`,
    which do not exist, "for `examples/ng_candidate_selection_probe.rs`", which does not exist.
  - `src/calling/allele_candidates/ssr.rs:1451-1454`: "Public because the measurement calls it
    … The run driver's dispatch will be the second caller." The example is gone; the driver
    (`src/run/callers.rs:1580`) is now the only caller.
  - Missing examples also cited at `likelihood/mod.rs:319` (`ng_minted_error_means.rs`) and
    `genotype_prior/seed_generic.rs:76` (`ng_seed_shape_weight_sweep.rs`).
  - Production paths that no longer exist (deleted 2026-09-13): `src/var_calling/…` and
    `src/ssr/cohort/…`, nine occurrences in this area (e.g. `allele_candidates/generic.rs:41`,
    `allele_candidates/mod.rs:1147`).
  - `likelihood/mod.rs:204` cites `tests::the_floors_still_equal_productions`; the test is
    `the_floor_still_equals_productions` (`:2932`).
- **Proposal** — fix each: delete the `hardy_weinberg` sentences, the nonexistent links and
  the "will be the second caller" sentence; for production paths and removed examples, say
  "the deleted production caller (commit `d9e7b076`)" or "a since-removed measurement program"
  so the history stays but no path points nowhere.
- **Downside** — none.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### R2. Public functions only tests use

- **Where**
  - `src/calling/likelihood/mod.rs:832` `fill_contaminant_allele_frequencies`, re-exported at
    `src/calling/mod.rs:89`. Production uses `fill_shared_contaminant_allele_frequencies`
    (`calling/mod.rs:2460`) plus `fill_own_batch_contaminant_row`
    (`inference/summarise_condition.rs:2336`); the all-in-one version is called only from
    `likelihood/mod.rs` tests (`:4041`, `:4466`, `:4482`).
  - `src/calling/likelihood/generic.rs:771` `genotype_log_likelihood_row`, re-exported at
    `src/calling/mod.rs:79` — see T1. It also shares its name with the production STR row
    `likelihood::ssr::genotype_log_likelihood_row`, so `calling::genotype_log_likelihood_row`
    is the generic one that production never calls, which is the opposite of what a reader
    would guess.
- **Proposal** — `#[cfg(test)]` both and drop the re-exports. The owner wants public
  interfaces minimal.
- **Downside** — the all-in-one contamination fill is a useful oracle for the split pair; keep
  it as a test helper rather than deleting it.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### R3. `count_deciding_partial_reads` copies a whole cohort observation to change one field

- **Where** — `src/calling/allele_candidates/deciding_partial_reads.rs:61-131`; the copy is
  built at `:124`.
- **What** — to add the deciding partial reads to each sample's rows, it clones the locus's
  allele table, its window coverage and its flank deletions, and builds a new
  `CohortObservation` whose only changed field is `per_sample`. The four selection passes then
  take the copy, and admission takes `merge_rows` separately. A reader has to work out that
  only `per_sample` differs.
- **Proposal** — return `Option<Vec<SampleSupport>>` and pass the per-sample rows to the
  selection passes as their own argument (they already take `merge_rows` that way).
- **Downside** — changes the signature of five `pub(super)` functions; the copy happens only
  at STR loci where some partial read singles out one allele, so the allocation saved is small.
  **Leave it** unless D1 is being done, which touches the same signatures.
- **Effect** — none (memory slightly lower). **Cost** — M. **Confidence** — medium.

### R4. `ssr_emission.rs` holds 2,730 lines of production code in one file

- **Where** — `src/calling/likelihood/ssr_emission.rs:1-2729` (tests from `:2730` to `:6530`).
- **What** — three separable topics: the single-segment model and its censored-read layout
  (`:430-1500`), the shared slip-placement and tract segmentation helpers (`:1500-1806`), and
  STR-bundle scoring (`bundle_emission` `:1807`, `bundle_censored_emission` `:2296`, and their
  helpers to `:2436`), followed by `CensoredLayout::lay_out` (`:2448-2729`).
- **Proposal** — split bundle scoring into `ssr_emission/bundle.rs` and the censored layout
  into `ssr_emission/censored.rs`, tests following their code.
- **Downside** — the owner's rule applies directly: the three share private helpers and
  scratch types (`MatchedPrefixLetters`, `Renderings`, `CensoredScratch`), which would become
  `pub(super)`, and a reader following one read's score would move between files. **Optional**;
  do it only if someone is about to work in the bundle code.
- **Effect** — none. **Cost** — M. **Confidence** — medium.

### R5. Doc comments carry history and measurement narrative at the expense of shape

- **Where** — throughout; examples: `select_generic` has an 80-line doc
  (`allele_candidates/generic.rs:18-80`); `AlleleSupport`'s doc (`run/cohort_merge/build.rs:1690-1749`)
  and `ReadGroupCalibration::charged_error`'s (`likelihood/mod.rs:491-545`) each include
  superseded positions ("Superseded in part on 2026-10-31", "which this sentence said until
  2026-08-25").
- **What** — the design principles ask type docs to lead with shape and add *why* only where
  it would surprise. Many docs here read as change logs: what a field used to mean, which
  review found what. The *why* is valuable; the superseded-position paragraphs make the
  current rule harder to find.
- **Proposal** — when a file is touched for another reason, move superseded positions and
  measurement narratives to the spec section the doc already cites, leaving one sentence and
  the citation. Do not do this as a sweep.
- **Downside** — this is the owner's established style and some of the history is what stops
  a reviewer re-proposing a rejected option. Ask the owner before any of it.
- **Effect** — none. **Cost** — L if swept, S per file. **Confidence** — low (a matter of
  taste).

---

## Looked at, left alone

- `GenotypePriorModel` trait (`genotype_prior/mod.rs:822`) with one production implementation
  and `&dyn` dispatch: test doubles (`CoefficientOnlyPrior`, `FixedLogPriors`) use the seam,
  and the dispatch is once per row.
- `fill_batch_allele_copies`' check that every batch holds a sample (`likelihood/mod.rs:747-755`)
  is batches × samples per locus per pass; today every run has one batch
  (`parameter_estimation/joint/sequencing_batches.rs:12-16`), so it costs one scan of the
  samples. Revisit only if a batching flag ships.
- `differ_by_one_indel` (`likelihood/generic.rs:425`) accepts one contiguous indel of any
  length; that is the intended meaning ("one insertion or deletion event").
- `deciding_partial_reads` counts nothing when any allele carries a flank deletion
  (`:64-77`), so `partial_read_fits`' flank-deletion handling runs only in tests; documented,
  measured, intended.
- `normalise_batch_rows` clamps each frequency to `[1e-12, 1]` without renormalising, so a row
  can sum to one plus a few times `1e-12`. Negligible.
- The misread fold in STR admission (`allele_candidates/ssr.rs:1196-1227`) allocates two
  `Vec`s per pair compared, but only for alleles carrying a flank deletion, which are rare.
- `SelectionScratch` destructuring with named `_` fields (`allele_candidates/mod.rs:701-724`,
  `generic.rs:107-120`): deliberate, so a new buffer is a compile error; keep it.
- `ssr_production_differential.rs` (test-only, 683 lines) re-implements production's three
  STR selection rules to compare against a recorded tomato fixture. Still a working regression
  check on ng's ladder and fold; whether to keep a production-parity test now that production
  is gone is the owner's call, not a code-quality one.
- Asserts kept in release across this area (`summarise_alleles`, `leftover_of`,
  `charged_error`): deliberate, per spec §8 — each guards a silent wrong genotype.
