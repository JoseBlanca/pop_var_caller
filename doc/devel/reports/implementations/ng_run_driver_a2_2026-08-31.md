# ng direct mode, step A2 — the construction checks

**Date:** 2026-08-31. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step A2. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md) §6.2, §7.1a;
[`../../ng/spec/parameters_file.md`](../../ng/spec/parameters_file.md) §6.
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §3.4, §5.
**Modules:** `src/ng/run/callers.rs`, `src/ng/run/mod.rs`, `Cargo.toml`.

---

## What landed

Four refusals at construction, three of them before a single file is opened, because each
condemns the whole run and opening a thousand files first would only make the message slower.

- **A cohort of no samples** (`RunError::NoSamples`).
- **Parameters assembled for another cohort** (`ParametersAreForAnotherCohort`) — one inbreeding
  coefficient per sample of this run, one calibration per read group.
- **The file-descriptor headroom** (`NotEnoughFileDescriptors`) — counted over *files*, and the
  message carries the count, the limit and `ulimit -n`.
- **And after the files are open, each sample's `@SQ M5` digests against the reference's**
  (`SampleAlignedToAnotherReference`), with the outcome recorded as an `AssemblyCheckOutcome`
  because *no sample was aligned to a wrong assembly* and *no sample could be checked* are
  different facts.

## Two of the three planned checks were already answered elsewhere, and finding that out changed the step

Neither was "already built here": one is built a layer up and one is built a layer down, and in
both cases A2 still had work left — a different check from the one the plan named. Investigating
before coding is what found the difference.

### The sample-name match cannot be done here at all

The plan asked for "the parameters' sample list matched against the run's **by name** (never by
position)". That match is real and it is already enforced — one layer up.
`ParametersFile::to_run_parameters_for` calls `refuse_if_not_this_runs_inputs`
([`bindings.rs`](../../../../src/ng/calling/parameters_file/bindings.rs)), which compares the
file's sample list against `ReadGroups::read_groups_per_sample()` by name and in order and
refuses naming the position where they diverge.

The caller cannot repeat it. What the caller receives is the *assembled* `RunParameters`, and
that type carries no sample names: a vector of inbreeding coefficients in the run's order, a
vector of calibrations in read-group order, and nothing that says whose they are.

What is still possible, and what nothing else prevents, is a `RunParameters` assembled for one
cohort handed to a caller opened over another — `open` takes the parameters and the read-group
table as two separate arguments. **So A2 checks the arity**, and its message says where the name
match lives so the next reader does not go looking for it here.

### The reference check was mostly the open gate's, and building it showed exactly what is left

The owner had already replaced the plan's original wording (the analysed regions agreeing across
samples — a psp fact, and vacuous in direct mode where one segmentation serves every sample) with
the contig agreement. **That, too, turned out to be built**: the open gate in
`SampleReads::open` compares each file's `@SQ` list against the reference's contig table — names,
lengths, order, **and the `M5` digests whenever the reference carries them**
([`open_bam.rs`](../../../../src/ng/read/input/open_bam.rs)).

The first version of the wrong-assembly test proved it the hard way. It built a BAM claiming
`ffffffffffffffffffffffffffffffff` for every contig, opened it against a reference read from its
FASTA, and expected the new check to refuse it. The file never reached the new check: opening it
failed with `md5 disagreement at index 0 (contig 'chr1')`.

**So `check_assembly` covers one case, and it is the ordinary one: a reference read from a
`.fai`.** That path hands back the contig table at once and verifies the FASTA on a background
thread, which is what keeps startup cheap — and it means the files open against a reference
carrying no digests at all, where the gate has nothing to compare. The digests exist only once
the run joins that thread. Comparing them then is what A2 builds, and it is why the caller takes
the *verified* reference as an input of its own rather than reading the one the files were opened
against.

The tests now open against the `.fai` arm and hand the FASTA arm's checksums over separately.
**That is equivalent for the comparison, not identical to a real run**: a run reads one reference
twice — the immediate `.fai` view and the finished view of the same FASTA — where the fixtures
build two FASTAs that are byte-identical because both come from the same contig list.

## A new dependency, and why it is not new

`rustix`, for `getrlimit(RLIMIT_NOFILE)`. Nothing in the tree read a resource limit and neither
`rustix` nor `libc` was declared — but **both were already in `Cargo.lock` as transitive
dependencies**, so promoting one costs no compilation that was not already happening: the lock
file grows by a single line. The manifest already carries this promotion twice, for the arrow
crates and for noodles' BAM and CSI, with the reason written beside each. `rustix::process` is a
safe wrapper, so no `unsafe` reaches this crate.

## The refusal that no fixture could otherwise reach

The descriptor limit a real machine reports is far above anything a test can build a cohort
against, so a check that read it inline would have a branch nothing could exercise — and an
unreachable refusal is one nobody has ever read the message of. The syscall is therefore split
from the decision: `refuse_without_descriptor_headroom` reads the limit and
`refuse_if_the_files_outnumber_the_limit` decides, so a test passes `Some(4)` and sees the
refusal, `None` and sees it pass.

## Verification

| claim | how it is held |
|---|---|
| a run given no alignment files is refused, and refused *first* | `a_cohort_with_no_alignment_files_is_refused`, whose parameters are another cohort's — because assembling parameters for a cohort of none panics inside the pre-pass, which is the failure this refusal comes before |
| parameters for another cohort are refused with both counts | `parameters_assembled_for_another_cohort_are_refused_with_both_counts` |
| a cohort's own parameters are accepted | `a_cohorts_own_parameters_are_accepted` — the check must not refuse what it is meant to let through |
| the descriptor refusal names the count, the limit and what to do | `a_run_that_would_run_out_of_descriptors_is_refused_naming_both_numbers` |
| enough headroom, and a platform reporting no limit, are not refusals | `headroom_above_what_the_run_needs_is_no_refusal` |
| the count grows with files, not with samples | `the_descriptor_count_grows_with_files_not_with_samples` |
| a sample aligned to another build of the assembly is refused, naming the sample and the contig | `a_sample_aligned_to_another_assembly_is_refused_naming_the_sample`, opened against a `.fai` reference so the check under test is the one that runs |
| a sample carrying the reference's own digests passes, and the run says how much it compared | `a_sample_carrying_the_references_digests_passes_and_the_run_says_what_it_compared` |
| a reference with no digests reports that nothing could be compared, rather than passing silently | `a_reference_without_digests_reports_that_nothing_could_be_compared` |
| the cheap refusals run before any file is opened | `the_cheap_refusals_come_before_the_files_are_opened`, whose one file would itself fail to open — so the two outcomes differ and the test can tell which check ran |

## ⚑ One decision is waiting for the owner, and one number is theirs to settle

**The empty-cohort refusal was built on the plan's instruction to raise it, not on a ruling.**
The plan says: "Raised at Checkpoint A; strike it if you would rather the command line refuse an
empty list." It is built and tested here so the choice is concrete, but it is still a choice.

**The descriptor arithmetic follows spec §7.1a and not this reader.** The spec says "a CRAM and
its index are two descriptors each"; read as built, an open `AlignmentFile` holds no descriptor
at all, its index is parsed into memory at open, and the descriptor belongs to a cursor — of
which direct mode holds one per file. So two per file over-estimates by about twofold today,
which refuses a run that would have fitted rather than letting one die at `EMFILE`: the safe
direction, and still the wrong number. It may stop being safe once several callers are in flight
(spec §11, question 2), which nobody has counted. The constant is left at the spec's value with
its documentation corrected to say what it is.

## What the reviews changed

Three reviewers read the step in separate worktrees. **The design was found faithful** — no
reserved decision taken, no B1 work leaked in, and all three of the rulings' premises verified
independently against the code. What they found was in the messages and in the reporting type,
which is most of what this step is.

- **The check gated on the wrong field.** It asked whether the reference had a *whole-reference*
  digest while what it compares is the *per-contig* checksums. Correct today only because the two
  travel together.
- **`Compared { files, contig_digests }` could report that nothing was compared** — files above
  zero and checksums compared at zero, which is a check that looked at nothing printed as a check
  that passed: exactly the substitution the type exists to prevent. And the count was summed
  across files with no denominator, discarding the one `AssemblyCheck` already offers. It is now
  `EverySampleMatchedTheReference { alignment_files, checksums_compared, checksums_possible }`,
  with `NothingCouldBeChecked { because }` naming which side had none — the reference or the
  files — since that is the side an operator can change.
- **The parameters refusal printed "1 samples"** at the single-sample end of the committed range,
  and carried a free-form string where the neighbouring error in the parameters module carries
  typed fields. It now reads "the number of samples is 2 in the parameters and 1 in this run".
- **The assembly refusal never named which reference the run was calling against** — two inputs
  in play, one of them wrong — and its wrapper repeated its own cause's sentence. It names the
  reference now, and a test asserts the diagnosis appears exactly once in the rendered chain.
- **The descriptor message quoted a total nobody could reproduce**: 34 where a reader computing
  "one file at two each" gets 2, because the 32-descriptor allowance was invisible. It now shows
  the arithmetic it is asking the operator to act on.
- **Three assertions could not fail for the reason they stated.** One checked that the message
  contained `'4'` — satisfied by the `4` in `34`, so it would have passed with the limit deleted
  from the message. All three now assert rendered messages.
- **The manifest cited the wrong precedent.** The arrow and noodles entries are about types this
  crate already named through another crate's re-export; nothing re-exported `rustix`. The
  fitting precedent is `serde_json`'s — already in the lock file, so it costs no new build — and
  even that needed narrowing: `tempfile` enables only `rustix`'s `fs` feature, so `process` does
  compile a module the build did not have.

**One thing the reviews retire rather than change.** `arch/run_streaming.md` §8 asked whether an
alignment header naming contigs the run does not analyse is a refusal or is ignored. It is a
refusal, decided by the open gate before A2 runs — `ContigList::first_disagreement` refuses on
any length difference — which is also what makes `check_assembly`'s parity assertion safe.

**And one methodological failure of my own, recorded because it nearly cost the step.** A
text-anchored bulk edit matched its start marker and not its end, and silently deleted 900 lines
of `callers.rs` including both test modules. The review worktrees held a byte-identical copy, so
the restore was exact and verified with `diff`, and the rewrite used line numbers instead. It is
the second pattern-matching edit in one session to do the wrong thing quietly — the first was a
mutation that reported "survived" when it had never applied — and both were caught by a check
made after the fact rather than by the edit itself.

## What the correctness review changed — and the gap it found that none of us saw

**Seven of twelve mutations survived, all in one blind spot.** Every fixture the first draft had
was one sample, one file, one read group — a shape in which *samples*, *files* and *read groups*
are the same number, so both arity checks fail together and neither is pinned; in which reversing
a two-element pairing is the identity, so the zip that decides **which sample gets blamed** for a
wrong assembly could be reversed unnoticed; and in which a file counter pinned to `1` reads as
correct. New fixtures separate all three: two samples with one read group each against parameters
for one sample with two (and the mirror image), one sample across two files, two samples sharing
one file, and a three-sample cohort with the wrong assembly on the *middle* one.

**And it found a gap in the step's own remit.** The plan's A2 asks for "each sample's alignment
header agreeing with **the segmentation's reference**". The step checked the samples against the
reference and never checked the segmentation against it — though both values were in hand.

The layer that should hold that check is blind on exactly the path A2 exists for.
`RepeatCatalog::open_checking_against_reference` guards both its per-contig and its
whole-reference digest comparisons behind the reference *having* a digest
([`reader.rs`](../../../../src/ng/repeat_catalog/reader.rs)), and a reference read from a `.fai`
has none. So on the ordinary path a catalog is admitted on contig names, lengths and order alone,
and nothing re-runs the comparison once the background FASTA read has finished.

**What that lets through is silent and genome-wide.** A repeat catalog built on one build of an
assembly, with a reference and alignment files on another build carrying the same contig names and
lengths: every sample's checksums agree with the run's reference, so the sample check passes; the
catalog's contigs agree on name and length, so its open passes; and every repeat tract's
coordinates are then applied at the wrong positions across the whole genome. `RunError::
CatalogIsForAnotherReference` closes it, comparing the catalog's own non-optional whole-reference
checksum against the reference's — two values already in hand at construction.

**A second refusal came with it.** Nothing tied the reference the samples were *checked* against
to the one their files were *opened* against, and the comparison downstream walks the two contig
lists in step — asserting that parity in debug and zipping in release. Two different genomes would
have paired a file's chromosome against something else's, blaming the wrong contig or missing a
real mismatch past the end of the shorter list. `ReferenceCheckedAgainstAnotherGenome` refuses it
on names and lengths before any checksum is read.

**One survivor is left, and it is untestable rather than untested.** Reading the *hard*
`RLIMIT_NOFILE` instead of the soft one survives, because the choice lives in the one-line wrapper
around the syscall and no test can control a process's real limit without `setrlimit`, which is
process-wide and would race a parallel test runner. The decision is documented at the call site
instead. Recorded rather than papered over.
