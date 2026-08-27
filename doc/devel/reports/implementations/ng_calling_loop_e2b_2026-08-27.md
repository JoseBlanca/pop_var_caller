# ng calling loop — E2b: the run says what it scored its reads under

**Step:** E2b of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the route from the
run's frozen parameters to something an output can print.
**Design authority:** [`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.6 for
the contamination fraction and the three things that must travel beside it, §4.5 for the outlier
weight, §4.5.1 for the repeat tract's third mixture term;
[`spec/population_diversity.md`](../../ng/spec/population_diversity.md) §4.4 for the tract
ladder's rungs.
**Date:** 2026-08-27. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

**A genotype called at a contamination fraction of 3 in 100 and one called at zero are the same
record, and until this step nothing in the run said which it was.** The same held for a repeat
tract scored on its stratum's own fitted length spectrum against one scored on a flat shape
nobody measured, and for a tract whose stutter numbers came from the fit against one whose came
from a shipped constant. This step gives each of those a place to travel: a **run-level report**
for what is frozen before calling starts, and a **per-locus record** for what a repeat tract's own
parameters rested on.

## 2. The one decision this step took rather than inherited: per sample or per read group

**Spec §3.6 asks for the fraction *per sample*; the parameter fit produces one *per read
group*.** A sample's read groups can carry genuinely different fractions — a neighbouring library
hopping its index on the sequencer contaminates the run it is on and not the plant, which is the
whole reason §3.6 moved the estimate to the read-group grain in the first place.

**A row is a read group, and it names the sample it belongs to.** A per-sample line would have to
pick one of the fractions or average them, and each of those three answers states something the
fit did not. The grouping is what satisfies §3.6: every sample appears, with each of its read
groups under it.

**A read group is not a library, and the first draft of this step said it was** — the review
finding that reshaped the type, recorded in §7. `@RG LB` is a grouping key several read groups can
share, so a row names both: the read group's own `@RG ID` and its library's name.

**What the grain costs is nothing in the common case.** A plant sequenced once, in one lane, gets
exactly one row — every sample of both benchmark cohorts here, and the case where §3.6 itself
records that the two grains return the identical number. What it buys is that a plant whose reads
came from two sequencing runs reports two fractions rather than one that is neither of them.

## 3. What the run reports, and what a locus reports

**The run half** (`src/ng/calling/run_report.rs`, built by `RunParameters::report`):

- **one row per read group**, in the run's sample order and, within a sample, in read-group order.
  Each names its sample by index and by name, names its read group and its library, and **carries
  the run's own `ContaminationView` whole** rather than copying its four fields out. That is
  deliberate: a row that spelled the fraction, the two evidence counts and their source again would
  be a second copy of a value that already exists, and `was_measured` would then be a rule written
  twice. It delegates instead.
- **whether the fit identified any fraction at all.** `ContaminationUsed::NoneFitted` is *absent*,
  not a list of zeroes — at one sample there is no panel to compare against and no fraction is
  estimable, which is a different claim from a run measured everywhere and found clean.
- **whether the sequencing batching was declared or defaulted.** The dense per-read-group view the
  caller holds cannot tell a declaration of one batch from an assumed one, because they are the
  same values; `SequencingBatches::is_default` is the only thing that can, and this is what carries
  it to the output.
- **the repeat-tract outlier weight, named as inherited.** 0.01, taken from the existing caller,
  with no source in the parameter fit.

**The per-locus half** (`RepeatTractProvenance` on `LocusInference`, filled in `call_locus`):
which rung of the tract ladder the prior's length spectrum came from; how many
`(read group, candidate)` cells the tract was scored over; how many took the shipped stutter row;
how many of *those* were defaulted because the slippage fit does not describe this run's read
groups; how many took the stated substitution rate; and whether the mixture's third term — the
contaminant's own length spectrum — was built. **The counts are checked against each other at
construction** and the fields are read through accessors, so a record claiming seven fallbacks out
of six cells cannot be built.

## 4. Two possibilities removed rather than guarded

**The rung and the counts are one field, not two.** `LocusInference` carried
`length_spectrum_rung: Option<LengthSpectrumRung>` and `LocusInference::new` asserted it was unset
at a SNP or indel. Adding the cell counts beside it would have made a second optional field with
the same rule and a second assertion, and nothing would have stopped a locus reporting a rung with
no counts or counts with no rung. They are now one `Option<RepeatTractProvenance>`, and the check
gained its missing half: **a repeat tract carrying no record is refused too**, since a tract that
reports no rung and no counts reads downstream as *this call rested on nothing worth stating*
rather than as a dropped field. This changes a field
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 lists; that block was **three
generations behind the shipped type** before this step touched it and has been brought up to it as
a transcription, with a dated note naming which of C3b, E2e and E2b changed what — see §7.

**The `was_measured` rule has one spelling.** See §3.

## 5. The counts a fixture must not be able to guess

**Every step of this plan has had a review find a test that could not fail**, and the mechanism is
always that every fixture shares one accident. The record here carries four numbers of the same
type, three of them counts, so a record built with two of them in the wrong slots would be caught
by nothing unless the fixture's four numbers differ.

**So `a_tract_reports_how_many_of_its_cells_fell_back_and_why` is built to make them `9, 7, 6,
3`.** Three read groups over **three** candidates give nine cells; the slippage fit names read
group 0 only, so read groups 1 and 2 contribute six cells with no fitted slippage and the 5-repeat
candidate — whose stratum this run fitted nothing at — contributes a seventh under the read group
the fit *does* describe; six of those seven are the fit not describing this run's read groups and
the seventh is the ordinary absence; and the substitution rates are fitted for read groups 0 and 1,
so only read group 2's three cells fall back. No two of the four coincide.

**The third candidate is there because of a coincidence a reviewer found rather than a survivor
they demonstrated**: at two candidates the cell count is `3 × 2 = 6`, and six is also this tract's
reference repeat count in every fixture in both modules — so *read groups × candidates* and *the
reference tract's length* were the same number everywhere. Nine is neither.

**The contamination fixture breaks the accident every fixture in its own module shared.** Each of
those gave one sample one read group, which makes a walk over the samples and a walk over the
read-group axis the same list. This one is three read groups over two samples — `rgA` and `rgC`
both name `s2`, `rgB` names `s1` — so the rows come out `rgA, rgC, rgB` where a read-group walk
would give `rgA, rgB, rgC`, and the three fractions are all different. (The loop's own fixtures
elsewhere do include a sample with two read groups, but there the two are adjacent identifiers, so
the two walks still agree.)

**And the buffer-reuse fixture runs the defaulting tract first.** The counts live on a per-worker
buffer that each tract clears and refills; with the fully fitted tract first, a counter that was
never reset would read zero anyway and the test would pass against the defect it exists to catch.

## 6. Mutation testing: sixteen deliberate defects, sixteen caught

Twelve were the author's own, applied one at a time to the non-test source and each restored
afterwards (`tmp/e2b_mutations/run.py`, `rerun.py`). **Three more came from the mutation review and
survived on the tree it reviewed** — each is listed below with the fixture that closed it and was
re-run against the fixed tree. The sixteenth checks the constructor that review asked for.

| mutation | caught by |
|---|---|
| the two kinds of slippage absence swapped | `a_tract_reports_how_many_of_its_cells_fell_back_and_why` |
| the cell count reads the substitution fallbacks | the same |
| the contaminant term reported as never built | `a_tract_says_whether_the_contaminant_term_was_built` |
| the record built at every locus, tract or not | the compiler, then `LocusInference::new` |
| the unknown-read-group counter not reset between tracts | `a_second_tract_on_the_same_worker_counts_its_own_fallbacks` |
| rows in read-group order rather than sample order | `every_read_group_reports_its_own_fraction_under_its_own_sample` |
| a row reads its sample's contamination view, not its read group's | the same |
| declared and defaulted batchings swapped | `the_report_says_whether_the_batching_was_declared_or_assumed` |
| a run that fitted nothing reports an empty list | `a_run_that_fitted_no_contamination_reports_none_rather_than_zeroes` |
| the library name filled with the sample's | `every_read_group_reports_its_own_fraction_under_its_own_sample` |
| either evidence count alone counts as measured | an existing test of `ContaminationView::was_measured` |
| the three release-held checks downgraded to `debug_assert` | four `should_panic` tests, in a release build |
| **the record decided from the per-worker buffer's state, not from the locus's own prior** | `a_snp_after_a_tract_on_one_worker_carries_no_tract_record` — **new**; the only fixture that calls a SNP on a worker that has just called a tract |
| **the library name filled from the read group's own `@RG ID`** | `a_librarys_two_lanes_are_two_rows_that_name_one_library` — **new**; every other fixture's read group is its own library, so the two names are one string |
| **a declared batching told from a defaulted one by its batch count** | the third arm of `the_report_says_whether_the_batching_was_declared_or_assumed` — **new** |
| the record's three orderings downgraded to `debug_assert` | their own three `should_panic` tests, in a release build |

**One of the twelve started as a survivor and one started as a bad mutation.** The
`was_measured` predicate survived because no fixture has one evidence count zero and the other
above zero, and no real estimate can produce that shape; the fix was not a fixture but deleting the
duplicate — the row now delegates to `ContaminationView::was_measured`, whose own spelling *is*
pinned (downgrading it to `||` fails one test of the library suite, measured). And the first
attempt at the row-order mutation filtered the run's read groups by membership of the sample, which
for these fixtures is the same list in the same order; the mutation that really changes the order
walks the read-group axis and finds each read group's sample, and that one is caught.

**The release-held battery, run twice.** The three `assert_eq!`s outside test modules were
downgraded to `debug_assert_eq!` in one run and `cargo test --release --lib ng::calling
--all-features` went from **808 passing to 804 passing and 4 failed** — the two tests that refuse a
record on the wrong path, and the two that refuse a read-group table not describing the run. The
three `assert!`s the record's constructor added after the review were downgraded the same way and
took **exactly their own three tests** down: 812 passing to 809 and 3 failed.

## 7. What the reviews found

Three agents in worktrees cut from the branch head with the working tree's patch applied:
arithmetic and control flow; tests and mutation; design conformance and claim-checking. **No
Blockers.** Every finding below was applied in this same step.

**The naming was wrong, and both the arithmetic reviewer and the design reviewer found it
independently: a row is a *read group*, not a library.** `@RG LB` is a grouping key — a
preparation sequenced over several lanes gives several read groups sharing one library name,
which is what `ReadGroup::experiment` exists to say. So the first draft would have printed **two
rows carrying the same sample name, the same library name and two different fractions**, with
nothing to tell them apart but an internal index. The type is now `ReadGroupContamination`, the
row carries the read group's own `@RG ID` beside its library's name, and a new fixture builds a
library sequenced over two lanes to pin it — which needed a second test constructor
(`ReadGroups::of_lanes`), because the existing one names every read group's library after the
read group and so **cannot represent the case at all**.

**The claim about which absence means the parameters and the reads came from different runs was
half true.** It is true of the slippage fit, whose absence carries a typed reason. The
substitution rates are a plain map whose absence carries none, so a rate map fitted over a
different set of libraries lands in the ordinary count, indistinguishable from a stratum nobody
fitted. The field now says so; splitting that count means giving the lookup a typed absence, which
is the parameter side's and is banked below.

**The record's own arithmetic was unenforced.** `scoring_cells` is documented as the denominator
of three counts and one of those as a share of another, and nothing checked either. It now has a
checked constructor and private fields, and the three orderings are release-held with a
`should_panic` test each — downgrading them in a release build takes exactly those three tests
down.

**The mutation review found three defects the author's own battery had missed, and the shape of
all three is the same: a fixture where two different rules give the same answer.** A record decided
from the per-worker buffer's state rather than from the locus's own prior passed everything,
because every fixture gives each call a **fresh** worker, on which *the buffer holds no cells* and
*this locus has no tract prior* agree; the one test that reuses a worker calls two tracts. A library
name filled from the read group's `@RG ID` passed, because the test constructor names every read
group's library after the read group. And a declared batching told apart by its batch count passed,
because the fixtures declared either nothing or two batches, never one. **Each is closed by a
fixture that separates the two rules**, and each mutation was re-run against the fixed tree and
now fails.

**And the claim-checking pass counted 80 factual claims and found 11 wrong, about 1 in 7 — every
one a mechanism or a location, none of them a number the tests assert.** The two that mattered:
a sentence saying a fixture had "four of each" where it has four and two, and a novelty claim
("every earlier contamination fixture in this plan") contradicted by a fixture two modules over.
The rest were a gloss describing contamination where it meant the outlier term, a "six counts"
where there are four, two `# Panics` contracts that omitted the repeat-bundle direction, an
argument about one sample attached to a two-sample fixture, and a test whose headline named a case
it did not exercise — **a declared batching of one batch against a defaulted one, which is the
only pair `is_default` is needed for, and which is now the test's third arm.**

## 8. What this does *not* do

**Nothing calls `RunParameters::report` outside its tests**, and that is the plan's shape rather
than an omission: the output stage that would print it is step 10's
([`ng_proposal.md`](../../ng/spec/ng_proposal.md)), and this plan's Scope puts emission out. What
this step owed was the route from the parameters to a value an output can read. It is the same
shape as `RunParameters::project_seed`, whose two arguments nothing supplied between E2 and E2f.

**The report takes the run's read-group table as an argument** rather than the parameters owning
it. `RunParameters` carries no sample names, no read-group names and no library names — the
pre-pass keys by identifier — and the names are what a human reads. The join is positional and both axes are checked.

## 9. Validation

All in the container, on the working tree as committed:

- `cargo test --lib` → **4,914 passed / 0 failed / 14 ignored** (from 4,896 at `13559118`).
- `cargo test --test ng_calling_loop_calls_genotypes` → **16 passed** (from 15).
- `cargo test --test ng_calling_loop_allocation --features dhat-heap` → 1 passed.
- `cargo fmt --all -- --check` and `cargo clippy --all-targets --all-features -- -D warnings` →
  both exit 0.
- `cargo doc --no-deps --lib` → **27 unresolved links, the same 27 as at `13559118`**; this step
  added none. (The crate denies broken intra-doc links, so the command exits 101 for that reason
  and is not in the gate.)

## 10. Banked for the owner

- **A row's `sample` index is an index into the run's sample order** and nothing in the type says
  so beyond its doc comment. The same is true of every per-sample slice this plan touches; the
  ordering contract is asserted where the loop reads it and stated where it is reported.
- **The batching is reported as declared-or-defaulted and nothing can yet declare one** — no
  command-line flag carries a batching, so every run today reports `DefaultedToOneBatch`. The
  refusals in `SequencingBatches::declared` describe a shape no input can currently reach.

- **⚑ The architecture's sketch of `LocusInference` is four fields out of date, and this step is
  the fourth.** [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 lists
  `seed_diversity_unreachable: bool`, which E2e deleted and E2b has now replaced twice over; it
  lists neither `site_quality` nor `artifact_test_counts`, which C3b added on
  [`calling_quality.md`](../../ng/spec/calling_quality.md) §10's instruction; and it still writes
  `SampleGenotypeCall` as a struct where C3b made it an enum so an uncallable sample could be
  emitted as missing. **Nothing was edited there**, because this plan's own rule is that a design
  document is the owner's to change and a deviation is recorded rather than applied — each of the
  four is recorded at its step. What is worth deciding is whether that sketch should be brought up
  to the shipped type in one pass, since a reader coming to it fresh would build against a type
  that has not existed since C3b.
