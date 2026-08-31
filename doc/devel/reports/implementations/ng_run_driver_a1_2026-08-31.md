# ng direct mode, step A1 — the caller object, constructed but inert

**Date:** 2026-08-31. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step A1. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md) §5.1.
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §1, §3.4, §5
(amended the same day, before this step, and again after it to record what landed).
**Modules:** `src/ng/run/segments.rs` and `src/ng/run/callers.rs` (both new), `src/ng/run/mod.rs`.

---

## What landed

The object a direct-mode run **is**, holding what a run holds, and nothing that iterates:

- **`Segmentation`** — the run's segments in genome order, beside `SegmentationInputs`, the
  record of the catalog, the repeat-tract criteria and the analysed regions they were computed
  from. Built once and shared, because the segments are a function of the reference, the
  catalog, the criteria and the requested regions and of no sample's reads, so k samples walk
  one list rather than k lists that happen to agree.
- **`AlignedFilesVariantCaller`** — one open `SampleReads` per sample of the run, plus the
  shared read-only state: the read-group table, the reference, the read filters, the
  segmentation, the fitted parameters, the calling loop's configuration, the
  candidate-selection configuration and the merge's parameters.
- **`RunError`** — the two failures A1 can produce: a catalog whose segment stream fails, and
  a sample whose alignment files will not open.

**19 tests** (11 in `segments.rs`, 8 in `callers.rs`). `cargo fmt --check`,
`cargo clippy --all-targets --all-features -- -D warnings` and `cargo test --lib ng::run` are
green in the container; the module's suite reports **308 passing, 0 failed**.

## The one decision this step actually makes: where the run's sample order comes from

The run's sample order is **the read-group table's first-seen order**, and there is one place it
is defined.

This matters more than its size suggests, and the plan says why: step D2 exists because three
different sample numberings meet in the calling loop — the merge names samples by index in the
run's order, the parameters name them by name, and the scratch rows are the run's order with
the uncallable ones closed up. A mismatch between them produces wrong genotypes rather than a
crash. A run that minted its own sample order beside the read-group table's would be building
the accident D2 is meant to catch.

`AlignedFilesVariantCaller::open` therefore opens **one `SampleReads` per entry of
`ReadGroups::read_groups_per_sample`, in that order**, and `sample_names()` reports that order.
Nothing in this step sorts, deduplicates or re-derives it.

## Departures from the architecture sketch, and why

Six, none of them a design change. The architecture was amended earlier the same day; its §3.4
sketch is what these depart from, and it has since been updated to the landed shapes.

**1. There is no `SampleInput` type.** The sketch's constructor took `samples: &[SampleInput]`.
The read layer already answers that question: a run builds one run-wide read-group table from
its file paths, and `SampleReads::open_only_sample`'s doc comment states the rule for every tool
that is not single-sample — "a tool that handles a cohort must not use it: it calls
`build_read_groups` itself and opens one `SampleReads` per entry of
`ReadGroups::read_groups_per_sample`". Minting a per-sample input type beside that table would
have created the second sample ordering the section above says must not exist.

**2. The four file-opening arguments are grouped as `AlignmentInputs`.** The table, the
reference, the read filters and the build-a-missing-index flag are one question — what this
run's reads are — and travel together into every sample's open. Grouping them also keeps the
constructor at six arguments rather than nine.

**3. `new` is `open`.** It opens, validates and index-checks every file of every sample before a
read flows, which is what the name should say. The architecture's `PspVariantCaller` already
uses `open` for the same reason.

**4. No concurrency knob.** The sketch's constructor took `callers_in_flight`. The owner ruled
on 2026-08-31 that Milestones A to D build against the single-threaded merge and that what
Milestone E becomes is decided from D3's measurement, so a knob with no consumer would be a
number nothing reads.

**5. `MergeParameters` carries four of the merge's five run parameters in three fields**, because
`MinAltReads` is itself a floor and a share. The one left out is how many building regions are
worked at once, which only means anything once the merge is threaded — the same ruling.

**6. `RunError` lands two variants where the architecture shows more.** The rest are A2's
refusals and psp mode's; a variant nothing constructs is a message nobody has written a test for.

## Three things added that the sketch did not name

**`Segmentation::analysed_regions()`, returning a slice.** The merge advances over the analysed
regions and takes them as `&[GenomeRegion]` (`merge_cohort_through_cache`); a sample's walk
advances over the segments. Both come out of the one object, built together in `build`, so the
two halves of a run cannot end up over different ground.

**`RunError::Catalog` carries the catalog's path.** Four of the catalog's own failures — a digest
mismatch, over-permissive criteria, differing scan weights, a differing tool version — describe
the file without naming it, so a person with several catalogs on disk could not tell which one
spoke. `Segmentation::build` takes the path for this and no other reason.

**Hand-written `Debug` for both objects.** A derived one would print every open file, every
segment and every fitted number — megabytes for a real cohort — in a line someone is reading to
find out which run this is. Both print names and sizes, and both have a test asserting the
contents do not appear.

## What is deliberately not here

- **The construction checks** — the parameters' sample list against the run's, each sample's
  alignment header against the segmentation's reference, the file-descriptor headroom, and (from
  this step's review) a cohort of no samples. Those are A2.
- **Any iteration.** No merge is driven, no locus is called, no record is written.
- **Any re-validation of the segment stream.** That a segment never crosses a contig and is
  never cut is the typed-region generator's guarantee, consumed rather than re-checked
  (spec §4.3).

## What the review changed

Three reviewers read the step in separate worktrees: one for correctness with mutation testing,
one for fidelity to the spec and the amended architecture, one reading the error messages and
public names as the person who runs the tool. **The design was found faithful — no reserved
decision was taken here.** What they found was a suite that could not fail, and a number that was
wrong.

**Six deliberate defects survived the first draft's tests, and all six were the same shape: a
fixture that passed a default and asserted the default back.** Such a test cannot tell "the
object held what it was given" from "the object replaced it with the default". The one that
mattered most: dropping the read filters and substituting the shipped defaults produced no test
failure *and no compiler warning*, because the argument still had a second use. A run would then
have called every position under the shipped thresholds while a person read their own numbers off
the command line — wrong genotypes, no failure.

Every fixture now differs from its type's default, and the assertions are paired with an
`assert_ne!` against the default so a future edit cannot quietly make them tautological again.

**Two more holes of the same kind:**

- `first_difference`'s field order was untested, because no fixture differed in two fields at
  once — the only situation in which the order is observable. Two tests now differ in two and
  three fields, and the order is documented as the order a person should fix them in: the catalog
  carries the reference's identity, so under a different catalog the other two comparisons are
  about different genomes.
- The sample-order test claimed "not the order the paths were given" and could not show it: with
  one sample per file, first-seen order *is* path order. There is now a fixture of one BAM
  declaring two read groups naming two samples, so the order is decided inside a file and neither
  path order nor alphabetical order can explain it.

**The wrong number.** The doc comment read "one open `SampleReads` per sample — 11 to 15 MiB of
live heap each". Spec §5.1 measures that **per open alignment file** (slope 12.0 MiB a file), so a
sample sequenced across four lanes costs four times what the comment claimed, and the 0.9 GB at 63
samples it quoted holds only at one file a sample. Corrected, with the multiplier named.

**Also corrected:** a summary line that said "variants in, genome order out" when variants are the
output; a paragraph whose stated reason for the shared reference said the opposite of what it
meant; `MergeParameters`' doc miscounting what it leaves out and calling a memory knob a thread
count; the claim that the catalog's error "names the file and the row", which names neither; and
three refusal strings that rendered as "written under a different routing criteria" — a word
naming nothing a person passes.

**One accessor was missing and one field was not held.** The read-group table was borrowed at
construction and dropped, and the parameters file the run must write beside its output is keyed by
read-group identifier with no table of its own — so writing it would have needed a second table,
which is the second-numbering accident again. The caller now holds it. `samples()` was added
because the next two steps both walk every sample.

## Verification

| claim | how it is held |
|---|---|
| a cohort of three opens, and a caller that sorted its samples would fail | `a_cohort_of_three_opens_and_names_its_samples`, whose sample names are not in alphabetical order |
| the run's sample order is the read-group table's, decided inside a file rather than by the path list | `the_sample_order_is_decided_by_the_read_group_table_not_by_the_paths` |
| two files of one sample are one sample | `two_files_of_one_sample_are_one_sample` |
| an unopenable file refuses at construction, and the person is told the sample, the file and that the index is missing | `a_sample_whose_index_is_missing_is_refused_naming_the_sample_and_the_file`, which asserts the whole rendered cause chain, not the wrapper alone |
| every setting is held as given, and none of them is its default | `every_setting_comes_back_as_it_was_given` |
| the shipped merge defaults are each type's own, not numbers restated here | `the_default_merge_parameters_are_each_types_default` |
| the ground and the read-group table come back | `the_ground_and_the_read_group_table_come_back` |
| neither object's debug rendering prints its contents | `the_debug_rendering_is_names_and_sizes`, `the_debug_rendering_of_a_segmentation_is_sizes` |
| a segmentation records the inputs it was given, none of them a default | `build_records_the_inputs_it_was_given` |
| a failing catalog stream fails the build, naming the catalog and carrying its reason | `a_failing_catalog_stream_fails_the_build_naming_the_file` |
| each differing field is named on its own | four tests, one per field plus the identical case |
| and when two differ at once, the one a person should fix first is named | `the_catalog_is_named_before_the_criteria_and_the_regions`, `the_criteria_are_named_before_the_regions` |

## Mutation testing, after the fixes

Every defect the correctness review found surviving was re-injected against the fixed suite, in
the review worktree, one full test compile each. **All seven are caught.**

| deliberate defect | before | after |
|---|---|---|
| `first_difference` checks its three fields in reverse order | survived | caught |
| the constructor stores `MergeParameters::DEFAULT` instead of its argument | survived | caught |
| it stores `CandidateSelectionConfig::DEFAULT` instead of its argument | survived | caught |
| it stores `ReadFilterConfig::default()` instead of its argument (no compiler warning either) | survived | caught |
| `MergeParameters::DEFAULT` restates a number rather than delegating to the merge's own | survived | caught |
| `Segmentation::build` records `StrRepeatCriteria::default()` rather than what it was given | survived | caught |
| the catalog's path is replaced by a placeholder in the failure | — | caught |

**One methodological correction, recorded because it is the failure mode this kind of testing
has.** The last row first reported as *survived*, and it had not been tested at all: the
substitution pattern was written against the source before `cargo fmt` reindented the block, so
it matched nothing and the run measured the unmutated code. A mutation that does not apply
reports exactly like one the tests cannot catch. The script now diffs the file against its
backup and refuses to report a verdict when the two are identical.
