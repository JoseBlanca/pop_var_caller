# Code-quality map: `src/calling/inference/` and the top-level files of `src/calling/`

Reviewed 2026-11-26 on branch `code-quality`, by reading the code only (nothing built or run).
Files covered: `calling/mod.rs`, `evidence_shaping.rs`, `genotype_table.rs`,
`genotype_table_parity.rs`, `loop_parity.rs`, `quality_parity.rs`, `run_parameters.rs`, and
`calling/inference/{mod.rs, summarise_condition.rs, repeat_tract_parameters.rs, slippage_refit.rs}`.
Not covered: `parameters_file/`, `likelihood/`, `allele_candidates/` (other reviewers), and
`genotype_prior/` and `quality/`, which nobody assigned to this review.

## Summary

1. **The two-way dependency between `calling` and `parameter_estimation` comes from one file.**
   Every compile-time import from `parameter_estimation` into `calling` goes through
   `calling/likelihood/run_context.rs`: the table of misread rates after a long run of one base.
   That file depends only on `float` and `types`. Moving it to a top-level module breaks the
   cycle at no cost. The other ~200 references, from `calling` to `parameter_estimation`, point
   the right way: calling reads fitted parameters.
2. **Calling reads its input types from `crate::run`.** `CohortObservation` and the types it
   holds are defined in `run/cohort_merge/build.rs`, and 9 production files in `calling/`
   import them. Moving the type definitions (not the merge logic) below both modules fixes
   this. The work spans three review areas, so it needs one coordinated change.
3. **`summarise_condition.rs` (10,664 lines) should become a directory of four files** along
   seams already in the code. Production is 2,810 lines; tests are 7,850 lines, and the tests
   are where most of the reading cost is.
4. **`calling/mod.rs` (5,848 lines, 3,040 production) holds six unrelated types**, the worker
   scratch alone ~1,640 lines. Split it into files and keep the re-exports. Its module comment
   claims "each imports downward", which the file itself contradicts.
5. **Code for features that were never built.** Two of the four allele-discovery modes cannot
   run, a `loop` runs exactly once (and a lint suppression admits it), and comments describe a
   "final prune" that does not exist. Delete. Separately, the per-locus slippage re-fit (about
   600 production lines) was measured as a near-no-op and kept switched off. **Keeping or
   deleting it is the owner's call.**
6. **The three `*_parity.rs` files are frozen regression snapshots now, not comparisons.**
   Keep the assertions, but rewrite their ~250 lines of narrative about a caller that no
   longer exists. Consider replacing the 732 KB genotype-table fixture with an independent
   formula.
7. **No likely bugs found.** Section 1 lists the paths traced and two minor stale assertions.

---

## 1. Bugs

**No bug with a concrete wrong-output scenario was found.** These paths were traced end to end:

- the order of rows against samples (`claim_row_for` → `run_sample_of_each_row` →
  `summarise_final_pass`'s per-sample check, `summarise_condition.rs:1325-1378`);
- the convergence denominator at one sample and with uninformative samples
  (`run_frequency_loop`, `summarise_condition.rs:807-910`);
- the leave-one-out restore in the slippage attribution (`attribute_reads_under_the_posteriors`,
  `:2066`);
- the contaminant-table lifecycle (`mod.rs:2323-2530`);
- the choice of primary alternative and its tie-break (`:1126`);
- the genotype-table cache (`genotype_table.rs:174`);
- the rounding of fractional read shares on the STR path (`evidence_shaping.rs:612`), which is
  exact because the merge never composes STR reads across records (`build.rs:698-700`).

### 1.1 A stale invariant comment in `LocusInference::new` describes STR bundles as unreachable

- **Where** — `src/calling/mod.rs:3560-3565`.
- **What** — The comment says a repeat bundle "could not honestly build" a repeat-tract record
  and that "nothing can reach it: `tract_candidates` refuses a bundle". Bundles are now called:
  `SsrLocusShape::Bundle`, `bundle_candidates`, `gather_for_bundle`, and `call_locus` builds a
  record for a bundle from `BundlePrior::rung()` (`summarise_condition.rs:2569-2578`). The
  assertion still holds. Only its explanation is wrong, and a reader trusting it would believe
  bundles are not genotyped.
- **Proposal** — Rewrite the comment to say that a bundle carries a record whose rung comes from
  its most variable run.
- **Downside** — None.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### 1.2 Comments and assertion messages cite a "final prune" that does not exist

- **Where** — `src/calling/mod.rs:126-127, 3458, 3492, 3594, 3624`; `inference/mod.rs:676`;
  `run/callers.rs:1775, 1784`.
- **What** — No function prunes the candidate table after the loop. `grep -i 'fn prune'` finds
  nothing, and `CandidateAlleles` has `admit` but no removal. The reasons given for
  `CandidateAlleles` being owned, and two panic messages, rest on this feature. A developer
  chasing one of those panics would look for code that is not there.
- **Proposal** — Reword to the real reason: the loop takes ownership so it can hand the table
  back inside `LocusInference`. Drop "final prune" from the panic texts.
- **Downside** — If the prune is ever built, these sentences come back.
- **Effect** — none. **Cost** — S. **Confidence** — high.

---

## 2. Architecture

### 2.1 Break the `calling` ↔ `parameter_estimation` cycle by moving `run_context.rs` to the top level

- **Where** — `src/calling/likelihood/run_context.rs` (547 lines). It is imported by
  `src/parameter_estimation/misread_after_runs.rs:17-20` and
  `src/parameter_estimation/joint/census_file.rs:57` (plus a test import at `:1985`).
- **What** — These are the only compile-time imports from `parameter_estimation` into
  `calling`. The other back-references (`joint/stratum_fits.rs:408, 416`,
  `joint/sequencing_batches.rs:190`, `joint/fit.rs:530`) are rustdoc links. The module holds
  one rule used by both sides, and its own doc says so (`run_context.rs:12-14`): which "entry"
  a base after a run belongs to.
  - The pre-pass counts by that rule: `contexts_along`, `ContextAfterRun`.
  - Calling charges by it: `letters_after_runs`, `a_run_may_reach_a_charged_length`,
    `ReadDirection`.
  - Both use the table types: `MisreadTable`, `RunEntry`, `ENTRY_COUNT`, `MisreadExcess`.

  The module imports only `crate::float` and `crate::types`, so it is a leaf.
- **Proposal** — Move it to `src/misread_after_runs_model.rs`, or a similar name beside
  `tandem_repeat.rs` and `genetics.rs`. Fix the 11 import sites. `calling/likelihood` and
  `parameter_estimation/misread_after_runs` then both import downward.
- **Downside** — Someone reading the STR likelihood finds the context rule one directory
  further away. Folding it into `parameter_estimation::misread_after_runs` would also break
  the cycle, matching how calling already imports `Slippage` from the fit. But a rule the
  scorer applies at every read then lives in the fitting code. The neutral module is the
  smaller compromise.
- **Effect** — none. **Cost** — S. **Confidence** — high.

**The other direction is not a problem.** The ~211 references from `calling` into
`parameter_estimation` are calling reading fitted parameter types: `Estimate`, `Provenance`,
`Slippage`, `StratumFits`, `StratumKey`, `SequencingBatches`, `ContaminationSource`. That is
the right direction. About 80 of those lines are in `calling/parameters_file/`, which
serialises the fitted parameters. Whether that sub-module belongs under `calling` at all is a
question for the `parameters_file` review.

### 2.2 Calling's input types live in `crate::run`

- **Where** — `src/calling/evidence_shaping.rs:77` (`CohortObservation`) and
  `src/calling/inference/mod.rs:67` (`MinAltObs`, `MinAltReadShare`, `MinAltReads`). The same
  pattern appears in 7 more production files in other areas:
  - `allele_candidates/{mod.rs:44-46, generic.rs:15, ssr.rs:46-47, deciding_partial_reads.rs:47, ssr_production_differential.rs:58}`;
  - `likelihood/{mod.rs:142, generic.rs:13}`.
- **What** — The cohort merge's output is defined in `src/run/cohort_merge/build.rs:1289-1800`:
  `CohortObservation`, `SampleSupport`, `SupportedAllele`, `PartialObservation`,
  `AlleleSupport`, `LocusBorder`. The support thresholds are in `run/cohort_merge/mod.rs`
  (`MinAltObs` at `:574`, `MinAltReadShare` at `:627`, `MinAltReads` at `:725`). These types
  are the interface between the merge and calling. The types themselves depend only on
  `locus_generation`, `types` and `window_coverage`. Only their constructors
  (`CohortObservation::over`, `over_from_cursors`) need merge internals.
- **Proposal** — Move the type definitions and their small accessors (`flank_deletion_of`,
  `pooled_support_for`) to a module both sides sit above, for example
  `src/locus_generation/cohort_observation.rs` or a top-level `src/cohort_observation.rs`.
  Keep the constructors as `impl CohortObservation` blocks in `run/cohort_merge/build.rs`
  (Rust allows inherent impls anywhere in the crate). Move the three `MinAlt*` types alongside.
  Re-export from `run::cohort_merge` so `run`'s own code does not change. **This touches three
  review areas, so do it as one change.**
- **Downside** — A reader looking for how a `CohortObservation` is built must find the `impl`
  block in another module. Name the constructor's location in the type's doc comment.
- **Effect** — none. **Cost** — M. **Confidence** — high.

### 2.3 `calling/mod.rs` describes a dependency rule it does not keep

- **Where** — `src/calling/mod.rs:5-11, 18-19`. Counter-examples: `:96-98` (imports from
  `genotype_prior`, `likelihood::ssr`, `quality`), `:109` (`likelihood::run_context`), `:1617`
  and `:2559, 2572, 2878` (`inference::repeat_tract_parameters::TractScoringFits`).
- **What** — The module doc says the shared vocabulary "needs nothing from the sub-modules" and
  that every sub-module "imports downward". In practice `CallingScratch` holds buffers from
  every sub-module, including the inference module's `TractScoringFits`. Within one crate this
  is not a build problem, but the stated rule is false, and a reader who trusts it will reason
  wrongly about what may import what.
- **Proposal** — Fold into finding 4.1: once `mod.rs` is split, rewrite the comment to say what
  is true. The vocabulary types (`CandidateAlleles`, `LocusEvidence`, `LocusInference`) import
  downward. `CallingScratch` gathers every sub-module's buffers by design, because a worker
  owns one of it.
- **Downside** — None.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### 2.4 Leave `LocusGenotyper` as it is (one production implementation)

- **Where** — `src/calling/inference/mod.rs:678`. The only production implementation is
  `summarise_condition.rs:2445`; there is one test fake at `inference/mod.rs:1109`.
- **What** — The trait exists for two other scoring methods that were never built (the
  architecture's "arms B and C"). `run/callers.rs` and `run/psp_caller.rs` are generic over it
  in 9 places.
- **Proposal** — **Leave it.** It costs one generic parameter, compiles to direct calls, and
  gives `run/` a narrow surface to depend on. Deleting it would make `run/` name the concrete
  loop type with its two type parameters, which reads worse.
- **Effect** — none.

---

## 3. Duplication

### 3.1 The STR prior is looked up twice per locus

- **Where** — `summarise_condition.rs:2569` (bundle prior, in `call_locus`) against `:1946`
  (the same `bundle_prior_of` call, in `fill_tract_likelihood_table`). For segments, `:2601-2606`
  (reference repeat count and `ssr_allele_length_variation_at`) against `:1668-1685`
  (`TractEvidence::prior`, the same lookup) called at `:1933`.
- **What** — Each STR locus resolves its allele-length-variation prior twice, by two pieces of
  code that must agree. For a bundle that also means two `Vec` allocations twice. If the slippage
  re-fit is on, the second lookup repeats on every round.
- **Proposal** — Resolve the prior once in `call_locus`, as an enum of `TractPrior` and
  `BundlePrior`, and pass it to `fill_what_no_pass_recomputes_for_round` inside
  `SlippageRoundBuild` or a sibling argument. Use the same value for the seed and the rung.
- **Downside** — One more parameter through a three-deep call chain that already takes seven.
- **Effect** — none (the same values are computed). It is slightly faster on STR loci.
- **Cost** — S. **Confidence** — high.

### 3.2 Repeated test helpers

- **Where** — `Ploidy::try_new(2)` appears 83 times in the crate. Seventeen files define their
  own `fn diploid()`, among them `loop_parity.rs`, `quality_parity.rs`, `run_parameters.rs`,
  `summarise_condition.rs` and `repeat_tract_parameters.rs`. FNV-1a is written out 7 times:
  `quality_parity.rs:96`, `loop_parity.rs:181`, and 5 places outside this area.
  `non_zero_default` is duplicated at `inference/mod.rs:213` and `run/cohort_merge/mod.rs:932`.
- **Proposal** — Add `Ploidy::DIPLOID` (a `const`) in `types.rs` and delete the `diploid()`
  helpers. Put one `#[cfg(test)] fn fnv1a` in a crate-level test-support module. Leave
  `non_zero_default` alone: two `const fn`s of three lines each.
- **Downside** — A crate-wide sweep touches many files for little gain. Do it only alongside
  other edits to those files.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### 3.3 Leave `RunParameters` and `FrozenParameters` as a pair

`RunParameters` owns the parameters; `FrozenParameters` borrows them, built by
`RunParameters::view`. Both carry about ten matching accessors. This is an owner and a view,
not duplication: `FrozenParameters` is `Copy` and crosses the worker boundary, and both sets of
accessors have callers (each is used 9 to 34 times). See "Looked at, left alone".

---

## 4. Readability

### 4.1 Split `calling/mod.rs` into one file per type family

- **Where** — `src/calling/mod.rs`, 5,848 lines: production `1-3040` (plus `3076-3684`),
  tests `3684-5848`. The production code falls into six families:

  | lines | content |
  |---|---|
  | 112-386 | `CandidateAlleles`, `ExpectedAlleleCopies` |
  | 386-814 | `GenericLocusSample`, `LocusEvidence`, `SsrLocusShape` |
  | 816-1354 | `FrozenParameters` (~540) |
  | 1355-2967 | `CallingScratch` and its six buffer-borrow structs (~1,610) |
  | 2969-3329 | `SampleGenotypeCall`, `RepeatTractProvenance` |
  | 3330-3683 | `LocusInference` |

- **What** — A reader who opens `calling/mod.rs` to learn what the calling module is meets
  1,600 lines of scratch-buffer management. `prepare_for_locus` alone is 440 lines including
  comments (`:1668-2110`).
- **Proposal** — Create `calling/candidate_alleles.rs`, `locus_evidence.rs`,
  `frozen_parameters.rs`, `scratch.rs` and `locus_inference.rs`, and move each family's tests
  with it. `mod.rs` keeps the module doc and `pub use` lines, so every `crate::calling::X` path
  still works. Private helpers (`resize_and_fill`, `locus_kind_word`) go with their users.
- **Downside** — Fields that `CallingScratch` reaches across (`allele_candidates::SelectionScratch`,
  `quality::SiteQualityBuffers`) are already `pub(crate)` paths, so nothing new becomes
  visible. Type-level private access between the families: `locus_kind_word` is used by both
  `LocusEvidence` and `LocusInference`, so it becomes `pub(super)`.
- **Effect** — none. **Cost** — M. **Confidence** — high.

### 4.2 Split `summarise_condition.rs` into a directory along its existing seams

- **Where** — `src/calling/inference/summarise_condition.rs`, 10,664 lines: production `1-2810`,
  tests `2811-10664`, 133 tests. The production code has four parts that already read as
  separate stages:

  | proposed file | lines now | content |
  |---|---|---|
  | `mod.rs` | 1-25, 2428-2808 | module doc, `SummariseConditionLoop`, `call_locus` (~400) |
  | `frequency_loop.rs` | 66-915 | `PassPrior`, `score_one_sample` (one sample's genotype probabilities), `sum_cohort_expected_copies` (the cohort summary), the convergence test, `run_frequency_loop` (~850) |
  | `final_pass.rs` | 960-1497 | `PooledArtifactCounts`, `pool_allele_reads`, `pick_primary_alternative`, `summarise_final_pass`, `mint_genotype`, `variant_class_of` (~540) |
  | `likelihood_table.rs` | 916-958, 1500-2427 | `TableReassembly`, `is_callable`, the evidence accessors, `weakest_warrant_at_the_locus`, the three fills, `attribute_reads_under_the_posteriors`, `ContaminantFrequencies`, `assemble_genotype_likelihood_table` (~970) |

  The tests already have section banners that follow the same seams: frequency-loop tests up
  to `:4735`, convergence `:4736`, final pass `:5486`, the driver and contamination `:6885`,
  cost invariants `:10386`. Their shared fixtures are `:2845-3230` and are reused across all
  sections.
- **Is it worth it?** Yes for the tests, less clearly for production. 2,800 production lines in
  one file is large but coherent: one algorithm, read top to bottom. The 7,850 test lines are
  what make the file hard to navigate, and they divide cleanly. If only one split is done,
  move the tests out (`summarise_condition/tests/{fixtures, frequency_loop, final_pass, driver, contamination, cost}.rs`)
  and leave production in one file.
- **Proposal** — Turn the module into a directory as tabled, with a `tests/fixtures.rs` holding
  the fixtures from `:2845-3230` and `:5489-5620`. Functions that cross files become
  `pub(super)`.
- **Downside** — `call_locus` then calls into three sibling files, and someone following one
  locus jumps between them. The module doc in `mod.rs` should list the files in calling order.
  Splitting production gives a smaller gain than splitting the tests.
- **Effect** — none. **Cost** — M (tests), M (production). **Confidence** — medium on the
  production half, high on the tests.

### 4.3 Remove the allele-discovery modes that cannot run, and the loop that runs once

- **Where** — `src/calling/inference/mod.rs`:
  - `DiscoveryMode::AgainstFrozenFrequencies` and `AgainstFullConvergence` (`:344, 347`);
  - `DiscoveryConfig::max_rounds` (`:397`) and `DEFAULT_DISCOVERY_MAX_ROUNDS` (`:204`);
  - `CallingLoopConfigError::DiscoveryNotBuilt` (`:650`) and its branch in `validate`
    (`:553-559`);
  - their tests (`:804-811, 835-836, 941, 1067`).

  Also `summarise_condition.rs:2455-2463` (the `#[expect(clippy::never_loop)]`) and the outer
  `let outcome = loop { … break outcome; }` at `:2651` and `:2765`.
- **What** — Both modes are refused by `validate`, so no run can select them. `max_rounds`
  is read only by a test (`inference/mod.rs:757`). The outer loop in `call_locus` exists to
  match pseudocode for those modes. It always breaks on its first pass, and a lint
  suppression with a nine-line justification keeps the compiler quiet about it. Together this
  adds about 120 lines and one level of nesting to the most-read function in the module.
- **Proposal** — Reduce `DiscoveryMode` to `Off` and `BeforeTheLoop`, delete `max_rounds`, the
  error variant and its tests, and remove the outer loop. Keep one sentence in the module doc
  saying round-wrapped discovery was examined and refused, with the research reference already
  there (`tract_genotype_accuracy_2026-09-03.md` §6.5).
- **Downside** — Bringing the design back would mean re-adding the variants. The research
  document keeps the reasoning, so nothing is lost.
- **Effect** — none (both modes are unreachable from the command line or the
  `NG_TRACT_DISCOVERY` variable, `cli/calling_run.rs:349-356`, which sets only `BeforeTheLoop` or `Off` at `:352, 355`). **Cost** — S.
  **Confidence** — high.

### 4.4 Decision for the owner: keep or delete the per-locus slippage re-fit

- **Where** — all of `src/calling/inference/slippage_refit.rs` (309 production lines);
  `SlippageRefitConfig` and its constants (`inference/mod.rs:119-160, 226-320` (the struct at `:257`));
  `attribute_reads_under_the_posteriors` (`summarise_condition.rs:2066-2165`); the round driver
  in `call_locus` (`:2630-2768`); `RefitEmissionCache` plumbing in
  `fill_tract_likelihood_table`; `TractScoringFits::{apply_refitted_slippage,
  effective_slippage_of_cell, frozen_slippage_cells}` (`repeat_tract_parameters.rs:984-1063`);
  `slippage_refit_rounds` on the scratch; two environment variables
  (`cli/calling_run.rs:270-291`, constants at `:285, 291`).
- **What** — The re-fit is off by default (zero rounds) and can be switched on only through
  `NG_SLIPPAGE_REFIT_ROUNDS`. The tract-accuracy report measured it on the tomato tracts:
  **7 genotype changes in 6,993 tracts, net +1 tract**, at the designed settings
  (`doc/devel/ng/research/tract_accuracy_program_report.md:717-727`). The verdict there:
  "discard as a default … the machinery stays built … the env switch remaining
  experiment-only". It costs about 600 production lines plus their tests, and it is the main
  reason `call_locus` has a nested round loop.
- **Proposal** — **Recommend deleting it**: the measurement is done and the result was close to
  nothing. If the owner wants to keep it, keep it as is; it is well isolated.
- **Downside** — Deleting removes two environment variables that `doc/user/options.md:61`
  documents (as "for development, do not set"). That document would need a one-line edit.
  Re-measuring on a different cohort (high depth, other species) would mean rebuilding the
  feature. The measurement covers tomato at about three reads per position only.
- **Effect** — none on default output. It removes a development switch. **Cost** — M.
  **Confidence** — high on the facts; the decision belongs to the owner.

### 4.5 The three `*_parity.rs` files: keep the assertions, rewrite the narrative

- **Where** — `src/calling/loop_parity.rs` (968 lines, 9 tests), `quality_parity.rs` (881
  lines, 15 tests), `genotype_table_parity.rs` (309 lines, 4 tests) with
  `testdata/genotype_tables_production.txt` (732 KB).
- **What** — Each one compared this caller against the deleted one. Since the deletion steps,
  each compares against **numbers recorded from the deleted caller at commit `d9e7b076`**,
  keyed by a digest of the inputs. So they are regression snapshots against an independent
  implementation: for the frequency loop and the site quality, the only tests that pin the
  numbers to something this code did not compute itself. That is worth keeping. But:
  - their module docs describe the comparison "as it was built, with production running"
    (`loop_parity.rs:1-93`, `quality_parity.rs:1-75`), around 170 lines of prose about
    functions that no longer exist, with links into deleted files (`loop_parity.rs:36, 73`);
  - `loop_parity.rs:899-923` compares the shipped defaults (`1e-3`, `50`, `ALPHA_REF = 1.0`)
    against literals typed beside them. In effect it is a test that these defaults never
    change, and it says so only indirectly.
- **Proposal** —
  1. Rename to `loop_frozen_answers.rs`, `site_quality_frozen_answers.rs` and
     `genotype_table_frozen_answers.rs`, or keep the names, and replace each module doc with
     about 15 lines. Say what is pinned, where the frozen numbers came from (the commit), and
     what to do when a deliberate change breaks one: re-derive the expectation independently,
     never re-record.
  2. Rename the defaults test to say what it guards: the frozen answers were recorded under
     these defaults.
  3. Optional: replace the 732 KB genotype-table fixture with an independent computation in
     the test. The VCF genotype ordering has a closed form, and the log multinomial
     coefficient can be computed as an exact sum of `ln k`. That removes the largest file in
     `src/`.
- **Downside** — Step 3 loses a bit-for-bit pin on the coefficients against the old code.
  The output-identity check being built would catch any change that matters to output.
  Rewriting the prose loses the history of why each divergence exists. Keep the two
  divergence explanations (the inbreeding correction and the pass count off by one); they
  are still true and still asserted.
- **Effect** — none (test-only). **Cost** — S for 1-2, M for 3. **Confidence** — high on 1-2,
  medium on 3.

### 4.6 Small items

| What | Where | Proposal | Cost |
|---|---|---|---|
| Two test-only re-exports nobody uses, defended by a comment that is wrong: deleting `pub mod likelihood;` already breaks the non-test `pub use likelihood::…` lines at `:80-91`. | `calling/mod.rs:70-79` | Delete both lines and the comment. | S |
| Test accessor named like the scoring function it returns the result of. | `calling/mod.rs:3102` `SampleGenotypeCall::score_best_genotype` | Rename `genotype_quality`. | S |
| `generic_beside_reference` builds the value through `generic()` and then matches it with an `unreachable!` arm. | `calling/mod.rs:552-566` | Make `generic()` call `generic_beside_reference(…, ReferenceBesideLocus::NONE)` and move the assertion there. | S |
| A 17-line comment about a type deleted on 2026-08-27, carrying one live to-do. | `run_parameters.rs:74-89` | Move the to-do (the soft count of segregating positions, spec §6.2) to the plan that owns it, and delete the comment. | S |
| About 20 file citations into the deleted caller (`src/var_calling/…`, `src/ssr/cohort/…`), some as broken rustdoc links. | e.g. `inference/mod.rs:80, 154`; `slippage_refit.rs:15, 25, 45, 54, 58, 118, 140, 210`; `summarise_condition.rs:93, 594, 2032, 2722, 8152`; `genotype_table.rs:11, 36, 498, 635`; `calling/mod.rs:1431` | Replace with "the previous caller, at commit `d9e7b076`, `<file>`" so they stay checkable. | S |
| Test section banners use plan step codes (`C2`, `C3b`, `D1`, `D2`), and so does a lint justification ("D1 has in hand"). | `summarise_condition.rs:4736, 5486, 6885, 10386, 1255-1260` | Keep the descriptive half of each banner and drop the code. | S |
| The module doc lists allele discovery as defaulting to `DiscoveryMode::Off`, but the shipped default is `BeforeTheLoop`. | `inference/mod.rs:19-25` (the variant's `#[default]` at `:333` is `Off`, while `DiscoveryConfig::DEFAULT` at `:404` uses `BeforeTheLoop`) | After 4.3, drop `#[default]` from the enum so that `DiscoveryConfig::DEFAULT` is the one source of the default. | S |

### 4.7 Possible speed-up on the STR path (measure before acting)

- **Where** — `src/calling/evidence_shaping.rs:542-608` (`SsrEvidenceScratch::rebuild`),
  called for every STR locus from `run/callers.rs:1607, 1681`.
- **What** — For each sample and each `(allele, read group)` row, `rebuild` clones the allele's
  bases (`:567`, a `Box<[u8]>`, so one heap allocation), and for each partial read it clones the
  bases and the witnessed positions (`:588-590`). The outer buffers are reused; these inner
  allocations are not. At 1,000 samples with three rows each, that is about 3,000 allocations
  and frees per STR locus. The merge built `CohortObservation` from `SequenceObservation`s, and
  this converts it back.
- **Proposal** — Measure first. If it shows up in a profile, let `SsrSampleEvidence` read the
  merge's rows plus the shared allele table directly, which removes the round trip. That
  change sits in `likelihood/` (another review area).
- **Downside** — `SsrSampleEvidence` currently has one input shape, used by tests and by the
  pre-pass. A second shape adds an abstraction.
- **Effect** — speed only, possibly. **Cost** — M. **Confidence** — low that it matters.

---

## Looked at, left alone

- `GenotypeTable::build` per locus: cached per thread up to ploidy 8 and 16 alleles
  (`genotype_table.rs:54-63, 174`), so it is not a per-locus cost.
- The relative tolerance of `1e-12` in `evidence_shaping.rs:93`: the two sums are both
  same-sign sums of logs, so the measured disagreements (about 1e-15 relative) leave three
  orders of magnitude of margin.
- `whole_reads` rounding of forward and placed-left shares for STR rows
  (`evidence_shaping.rs:612`): exact on the STR path, since the merge composes reads across
  records only at SNP/indel loci (`run/cohort_merge/build.rs:698-700`).
- Partial STR reads get `num_fwd: 0` (`evidence_shaping.rs:596`): harmless, because partials
  are cached with one emission for both strands, and the row then takes the single-sum branch
  (`likelihood/ssr.rs:769-775`).
- `score_one_sample` overwriting a sample's copies in the final pass and in the slippage
  attribution: the final pass is one extra E-step against the converged summary, as designed,
  and the attribution restores the copies (`summarise_condition.rs:2093-2157`).
- `RunnableCallingLoopConfig` with `Deref`: a validated newtype whose only constructor is the
  check. Fine.
- `RunParameters` (owner) and `FrozenParameters` (borrowed `Copy` view): both accessor sets are
  used; this is not duplication.
- The `assert!`s held in release builds throughout `CallingScratch` and `score_one_sample`:
  they are O(alleles) or O(genotypes) per sample and guard joins that would otherwise give
  wrong genotypes silently. Keep them.
- `repeat_tract_parameters.rs` living under `inference/` rather than `likelihood/`: it
  assembles fitted parameters into what the STR row reads, and only the loop and the scratch
  use it. Either home is defensible, and moving it gains nothing.
- `evidence_shaping.rs`'s module doc (63 lines) is long, but each paragraph records a join the
  types do not enforce. Keep it.
