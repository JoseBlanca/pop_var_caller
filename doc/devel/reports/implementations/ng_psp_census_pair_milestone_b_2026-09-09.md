# The census lives inside the psp — Milestone B: one judgement, and a cohort that agrees on its settings

**Date:** 2026-09-09
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone B
**Spec:** [psp_census_pair.md](../../ng/spec/psp_census_pair.md) §4.1, §4.2, §6, §9;
[psp_file_format.md](../../ng/spec/psp_file_format.md) §3.4, §6.5
**Branch:** `census-vs-psp-perf`, on top of `b87c297d`
**Milestone A's report:** [ng_psp_census_pair_milestone_a_2026-09-09.md](ng_psp_census_pair_milestone_a_2026-09-09.md)

---

## The tree is not green, and every step is judged against this and not against green

Measured on the untouched tree at `b87c297d` — the head of Milestone A plus B0 — before the first
line of B1 was written. It is Milestone A's baseline re-run, and it agrees with what that report
recorded.

| gate | on the untouched tree at `b87c297d` |
|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | **6,682 lib tests pass**, 15 ignored; 20 of 21 targets green; the one failure is `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` (`tests/ng_calling_loop_calls_genotypes.rs:1241`) |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | **11 errors of 5 kinds in 6 files**: `src/bam/alignment_input.rs`, `src/ng/window_coverage/accumulator.rs` (5), `src/ng/run/cohort_merge/build.rs` (2), `src/ng/run/cohort_merge/serial.rs`, `src/ng/run/psp_writer_line.rs`, `src/ng/window_coverage/production_parity.rs` |
| `cargo check --all-targets --keep-going` | **4 examples** do not compile: `ng_call_cohort_end_to_end` (1 error), `ng_candidate_selection_probe` (4), `ng_cohort_merge_parallel_cost` (8), `ng_cohort_merge_real_cost` (11) |
| `cargo fmt --check` | **4 files**: `examples/ng_call_cohort_end_to_end.rs`, `examples/ng_call_from_psps_cost.rs`, `examples/ng_census_locus_spans.rs`, `src/ng/psp/block.rs` |

**The gate each step is compared against is a *set*, not a count.** The comparison lists clippy's
error kinds and their project-file locations, the targets `check` cannot compile, and the files
`fmt` rewrites, and diffs each list against the row above — because a check written from the
baseline's own failures cannot see a failure it has never seen before, which is how a step in
Milestone A reported "unchanged" over two new lints of a kind the filter did not name
(`reporting-in-chat`'s 2026-09-07 entry). The script that does it is `tmp/b/gate.sh`, which is
scratch and not committed.

---

## B1 — the verdict type

**Committed:** see `git log` for `feat(ng): B1`.

### What it does

`CensusVerdict` is what a run has to say about the census inside one psp: it carries one, in the
format this build reads and against the loci this run rebuilds (`Fresh`); it carries none, because
its trailer is empty (`NoCensus`); it carries one of another version (`AnotherFormat`); its trailer
holds bytes that are not a census (`NotACensus`); or it carries one written against a different
selection (`AnotherSelection`). Each renders as the last column of spec §4.3's one line a sample —
`<sample>  <psp path>  carries no census`. Nothing produces a verdict yet; B2 and B4 do.

`VERSION`, the census layout's version word, becomes `pub` so the type can name it.

### Two departures from the plan's words, both recorded

- **`NotACensus` is a fifth verdict where spec §4.2 lists three causes.** The read that reaches
  §4.2's second cause — the version word at the trailer's front — has to know the bytes are a
  census before it can believe the word, and bytes that are not are damage rather than an old
  format. Spec §4.2's own point is that an old census must not be reported *as* damage; that
  requires damage to remain sayable.
- **The plan says "No logic", and the type carries one constructor**,
  `CensusVerdict::of_a_version_word`. It exists because the review found the alternative unsafe:
  see below.

### What the review found, and what was done

One read-only agent over five grouped categories, forbidden to edit or run `cargo`, asked to name
mutations rather than run them; the report is
[psp_census_pair_milestone_b_2026-09-09.md](../reviews/psp_census_pair_milestone_b_2026-09-09.md).
Two findings changed the type's shape and four corrected doc claims that the code beside them did
not support.

- **`AnotherFormat` carried both version numbers as fields, and the pair could disagree with the
  build it came from.** `AnotherFormat { in_the_psp: 4, this_build_reads: 4 }` was one struct
  literal away from any caller, and it rendered as *carries a census built by **a newer** version
  of this program: it is version 4 and this build reads version 4*. The variant now carries the
  psp's version alone and the message reads this build's from the constant, so the contradictory
  value cannot be written.
- **The constructor returned `Fresh` for a version word that matched.** `Fresh` is all three of
  spec §4.2's causes and a version word is one of them, so a `regenerate-census` that skipped on
  that answer (spec §8) would skip psps it owes. It returns `Option<CensusVerdict>` now, and
  `None` — nothing to report — is what a matching version word gets.
- **`NoCensus` and `NotACensus` both described an empty trailer**, so B2 could have implemented
  either and both would have looked right against the docs. `NotACensus` now says *a trailer that
  holds bytes*, and the two partition the psps.
- Four doc claims corrected: the count of causes ("the three that are not `Fresh`" — there are
  four, in three cost classes); the cost of reaching `AnotherSelection`, which was written as
  *seconds a sample* when the reference read and the selection rebuild are **once a run** and only
  a digest comparison is per sample; the past tense of *reached the user as malformed*, which is
  what `decode_census` still does today; and the claim that a trailer which is not a census leaves
  *the psp's records whole*, which the open has not checked — what it checked is the footer and the
  index.
- `VERSION`'s own doc said **"2 since 2026-08-16"** while the constant is 4. Read out of the
  history rather than guessed (`git log` over the file, printing the constant at each commit
  that touched it): 1 at creation on 2026-08-14, 2 on 2026-08-16 with the wider depth code, 3 on
  2026-09-05 when the census began recording its read groups, 4 the same day when it began
  carrying what the base qualities claimed.

### Four mutations, all run

| mutation | before the fixes | after |
|---|---|---|
| the message prints this build's version where it was handed the psp's | **green, 4 passed** | 2 tests fail |
| the older/newer decision reads the constant rather than the field | did not compile — proves nothing | — |
| the older/newer decision is the wrong way round | — | 2 tests fail |
| the constructor calls this build's own version a cause to report | — | 1 test fails |

The first is the review's finding measured rather than argued: both tests set the *this build
reads* field to `VERSION`, so a message that ignored the field entirely passed all four. The fix
deletes the field, and the tests now write the found version as a literal and this build's from the
constant, so a message printing either number twice fails.

Each mutation was reverted and the file compared against its backup before the next; the restored
tree's four tests pass.

### What was measured

- **`cargo test --lib --bins --tests --all-features --no-fail-fast`: 6,686 lib tests pass** against
  the baseline's 6,682 — this step's four — with 20 of 21 targets green and the one pre-existing
  failure unchanged.
- **`cargo clippy --lib --bins --tests --all-features -- -D warnings`: the same 11 errors of the
  same 5 kinds in the same 6 files as the baseline.**
- **`cargo check --all-targets --keep-going`: the same 4 examples, with the same error counts.**
- **`cargo fmt --check`: the same 4 files.**

---

## B2 — the cheap half of the judgement

**Committed:** see `git log` for `feat(ng): B2`.

### What it does

`what_the_footer_and_the_trailers_head_say_about_a_census` takes an open psp and answers two of
spec §4.2's three causes for one short read or none. The footer says how long the trailer is and
was read when the file was opened, so an empty trailer is *no census* without touching the file.
Otherwise the first ten bytes of the trailer are read — a census's magic and the version word
behind it — and they are either not a census at all, a census of a version this build does not
read, or this build's own.

Two pieces underneath it. `PspReader::trailer_head(at_most)` gives the front of the writer's
closing payload, and `trailer()` is now that with no limit. `census_file::version_word_of(head)`
says which version a stretch of bytes claims, or that it is not a census — which is what lets an
old census be named as old rather than as damage, since `decode_census` refuses both with the same
`Malformed`.

**A cheap read is only cheap if something measures it**, so `PspReader` gains a per-thread count of
the trailer bytes it has handed out — the instrument the census's own reader already carries for
the same argument (`census_file::bytes_read`). See the mutation table: without it, the step's whole
property was unpinned.

### What the review found, and what was done

One read-only agent over five grouped categories; the report is
[psp_census_pair_milestone_b_2026-09-09.md](../reviews/psp_census_pair_milestone_b_2026-09-09.md).
It walked every input a psp can present — an empty trailer, one of 1 to 9 bytes, a wrong version, a
version above this build's, a whole census — and found no wrong verdict. What it found was one
missing measurement, two false doc claims, and a cast.

- **Nothing failed if the judgement read the whole trailer.** Measured: with
  `trailer_head(BYTES_THAT_NAME_THE_VERSION)` replaced by `trailer()`, all 64 tests passed. The
  property the step exists for — a thousand-psp cohort does not pull a thousand censuses into
  memory — was asserted in three doc comments and measured nowhere. The trailer-byte counter and
  two tests close it.
- **The `# Errors` section named a failure that cannot reach the function**: *a footer that points
  past its own end*. `PspReader::open` refuses any psp whose trailer does not end exactly where the
  footer begins, so such a file never becomes an argument. The doc now says what is left — a file
  truncated or replaced after it was opened, or an I/O fault — and says that `open` is what rules
  the other out.
- **The container module claimed to know what is in a trailer**, two lines under its own sentence
  saying it must not. `psp_file_format.md` §3.4 keeps the payload the writer's business precisely
  so that adding to it is not a container version bump. The size argument stays, as an example of
  what ng's walk puts there; the ownership claim is gone.
- **`trailer_bytes as usize` wraps on a 32-bit target**, turning a trailer wider than a `usize`
  into a few bytes reported as the whole of it. `usize::try_from(…).unwrap_or(usize::MAX)`.
- **The module summary said "the verdict and nothing that reaches it"**, which B2 falsifies.
- **The shared fixture went in the wrong door.** `census_file`'s `mod tests` had been made
  `pub(crate)` so another module's tests could reach one function; the project's own answer is a
  `tests_support` module beside it, which `psp::writer` has and which this step's own tests already
  import from. The fixtures moved there and `mod tests` is private again.
- Renamed: the judgement was `what_the_psps_head_says_about_its_census`, and *the psp's head* reads
  as the psp's **header** to anyone holding the format spec — which is the one part it does not
  touch. It reads the footer and the front of the trailer, both at the file's tail.
- Smaller: `trailer_head`'s parameter is `at_most` rather than `bytes`; the constant sits beside the
  magic and the version it is made of rather than 530 lines away; the fixture helper's `sample`
  argument says what it will be for in B4; and a census of a version *newer* than this build's now
  goes through the judgement as well as through the message.

### Six mutations, all run

| mutation | before the fixes | after |
|---|---|---|
| the judgement reads the whole trailer, not its front | **green, 64 passed** | 1 test fails |
| the head read ignores how long the trailer is | — | 13 tests fail |
| the empty trailer is not answered out of the footer | — | 2 tests fail |
| the two damage verdicts are swapped | — | 4 tests fail |
| this build's own version word is not read as fresh | — | 2 tests fail |
| the version word is always read as this build's | 1 test failed, and **not the one in the module that owns the function** | 2 tests fail |

The last row is why `census_file` now has its own assertion on a word that is not this build's: its
test wrote and read a census of the current version only, so a `version_word_of` that answered
`VERSION` whenever the magic matched passed it, and only the psp-side test caught it.

**One mutation stayed green and the claim was changed instead of the test.** Removing the
empty-payload short-circuit inside `PspReader::trailer_head` leaves every test passing, because a
zero-length read reads zero bytes and the counter counts bytes. What the code avoids there is a
seek and a syscall, not a byte, so the doc now says *no byte of trailer is read* — which is what
the test measures — rather than *no read at all*, which nothing here measures.

Each mutation was reverted from a backup and all three files compared against it before the next;
the restored tree's 66 tests pass.

### What was measured

- **`cargo test --lib --bins --tests --all-features --no-fail-fast`: 6,696 lib tests pass** against
  the baseline's 6,682 — B1's four and this step's ten — with 20 of 21 targets green and the one
  pre-existing failure unchanged.
- **`cargo clippy --lib --bins --tests --all-features -- -D warnings`: the same 11 errors of the
  same 5 kinds in the same 6 files as the baseline.** The first run of this step was **not**:
  it had a twelfth, `writing &PathBuf instead of &Path`, in a test helper of mine — caught because
  the gate diffs the list of error *kinds* rather than counting the ones the baseline already had.
- **`cargo check --all-targets --keep-going`: the same 4 examples, with the same error counts.**
- **`cargo fmt --check`: the same 4 files.** The first run of this step had a fifth, this step's own.

---

## B3 — the cohort agrees on its settings

**Committed:** see `git log` for `feat(ng): B3`.

### What it does

`OpenPspCohort::open` compared the analysed regions across a cohort's psps and nothing else. It now
compares all three of the settings a segmentation is a function of — the catalog, the repeat-tract
criteria and the ground — through `SegmentationInputs::first_difference`, and refuses a cohort that
disagrees, naming both samples and the field. The ground keeps its own refusal, because its fix is
its own: re-walk one of the two, or call each over the ground it has. The catalog and the criteria
get a new one, `RunError::CohortWalkedUnderDifferentSettings`.

Every command that opens a cohort this way gets it: `call-from-psps` and `generate-census` today,
`estimate-parameters` at plan step C2, which is when it stops opening its cohort the other way.

### What this buys, stated correctly — the review corrected me here

**It is not the first check of these fields, and my first draft said it was.**
`PspVariantCaller::open` has compared every psp's catalog and criteria against **the run's**
segmentation since psp mode was built, and `call-from-psps` builds that segmentation from its own
flags — so a calling run already refused a cohort like this. What B3 adds is a refusal that **names
the pair** rather than one sample and the run, and a refusal for the commands that never build a run
segmentation at all: `generate-census` today, and the parameters fit after C2, which is exactly the
command that would otherwise take the criteria from the first psp and learn twenty seconds later
from a digest that the others disagree.

**The spec's own premise sentence is wrong on the same point.** `psp_census_pair.md` §6 says
`first_difference` "is called only from its own tests"; it is called from `psp_caller.rs`'s caller
check, and was before this branch. The decision the sentence supports is unaffected — the check
still belongs in the shared opener — so the spec is left as it is and this is raised at the
checkpoint.

### One test had to change, and that is a real consequence

`a_psp_walked_under_another_catalog_is_refused_naming_the_field` built a cohort whose *second* psp
carried another catalog, let the opener pass it, and asserted the caller refused it naming that
second sample. The opener now refuses that cohort, so the test cannot reach the caller. It is
rewritten with **both** files carrying the other catalog: the cohort agrees with itself and
disagrees with the run, which is what the caller's check is for now that a cohort reaching it always
agrees internally. The property it used to carry — that a per-file check reaches every file — is
still pinned, by the contig-table test, whose subject is not part of what the opener forces equal.

### What the review found, and what was done

- **Two doc claims of the form "nothing checked this until now" were false**, for the reason above.
- **The new error's doc named commands that do not have the check**: `estimate-parameters`, which
  opens its cohort through `open_census_cohort` and compares the regions alone, and
  `regenerate-census`, which does not exist yet. It now says what is true today and names the step
  that completes it.
- **The `call-from-psps` test's account of its own fixture was wrong.** It said the second walk's
  purity floor types no tract; the fixture's reference is one base repeated, so every tract on it is
  perfectly pure and clears either floor. What differs is the criteria record in the header, which
  is what the opener compares — which is a fine basis for the test and is now what it says.
- **That test's assertions could not tell the two refusals apart.** Without the opener's check the
  run reaches the caller's, which names the same field; only the pair of sample names separates
  them. It asserts the whole sentence now.
- **The precedence between the two refusals was unpinned.** A cohort differing in both the ground
  and the catalog is refused about the catalog — a person cannot act on *different ground* between
  files that are not about the same assembly. A test says so.
- The field name the routing matches on is a constant, `SegmentationInputs::ANALYSED_REGIONS`, and
  its own test now reads it rather than repeating the literal.
- The architecture doc's list of refusal axes said eight and named them; it says nine.

### Six mutations, and one that did not compile

| mutation | outcome |
|---|---|
| the check reverted to the ground alone — what it did before this step | 4 tests fail |
| the criteria comparison dropped, the catalog kept | 2 tests fail |
| only the first psp examined | 5 tests fail, including the ground refusal that predates this step |
| the routing swapped, so the catalog is refused as a ground disagreement | 3 tests fail |
| the field-name constant given a different value | 1 test fails — its own, and no other, which is the point of the constant |
| the caller's per-file loop compares the *first* psp every time | **all 74 pass** |

The last is not a defect and was run to settle a claim: after this step the caller's segmentation
check cannot differ per file, because the opener has already forced the three fields equal across
the cohort. That is why the rewritten test above no longer claims to prove the caller's loop reaches
every file.

**The routing mutation the review named did not compile** — the two error variants have different
fields, so exchanging their bodies is not a one-line edit. A mutation that does not compile proves
nothing, so it was replaced by one that does: routing the catalog to the ground's refusal.

Each was reverted from a backup and both files compared against it before the next; the restored
tree's 53 tests pass.

### What was measured

- **`cargo test --lib --bins --tests --all-features --no-fail-fast`: 6,700 lib tests pass** against
  the baseline's 6,682 — B1's four, B2's ten and this step's four — with 20 of 21 targets green and
  the one pre-existing failure unchanged.
- **clippy: the same 11 errors of the same 5 kinds in the same 6 files as the baseline.**
- **`cargo check --all-targets --keep-going`: the same 4 examples.**
- **`cargo fmt --check`: the same 4 files** — the first run of this step had two more, both its own.

---

## B4 — a cohort judged whole

**Committed:** see `git log` for `feat(ng): B4`.

### What it does

`what_the_heads_say_about_every_census_in_a_cohort` takes an opened cohort and returns one row a
sample, in the order the paths were given: the individual out of the psp's header, the file the run
was given, and what B2's two cheap reads say about that psp's census. **Nothing stops it.** A stale
census does not, a psp that will not read does not, and that is the whole of the step: the cohort
opener it replaces returns at the first census it cannot check, so a cohort with three stale ones
reports one, the user regenerates that one, and the next run reports the second. Regenerating a
census is a quarter of an hour a sample (spec §4), so what the early return costs is that wait once
a stale sample, in series, to learn a job that fits in one message.

A psp that cannot be *read* is not a psp whose census is stale. Its row carries the read failure
instead of a verdict — `Result<CensusVerdict, PspReadError>` — because regenerating the census of a
truncated file would not fix it, and a report that sent the user at it would be pointing at the
wrong fault. **The step's trap was writing that loop with `?`**, which reproduces for read failures
exactly the early return §4.1 exists to remove; B2's review named it in advance.

Nothing calls this yet. `estimate-parameters` takes its verdicts at plan step C3 and
`regenerate-census` at D2.

### Two choices the plan left open, both small

- **The row is a type and not a tuple.** `JudgedPsp` carries the sample, the path and the verdict,
  because spec §4.3's report is a line a sample with all three on it and because neither of the
  first two names the other: a psp's header carries the individual, its filename is whatever the
  walk was told to write, and the two agree only by convention.
- **The loop reads the cohort through a new accessor.** `OpenPspCohort::each_psp_with_its_path`
  hands out each open reader beside its path — the readers are private and a trailer read needs
  them mutably. The judgement itself stays in `census_freshness.rs` beside the verdict it produces,
  rather than moving into the module that owns the cohort.

### What the review found, and what was done

**The one that mattered: a property four tests claimed and none could see.** Every fixture psp was
written to `<sample>.psp` with its header naming that same sample, so *the row names the right
individual* and *the row names the file's stem* were the same assertion. A pass that read the
individual off the path — which is wrong on any cohort whose files are named by run accession or by
lane — would have passed all of them. **Fixed in the fixture**: a cohort's psps are now written to
`psp-0.psp`, `psp-1.psp`, … while their headers name `delta`, `alpha`, `echo`, `bravo`, `charlie`,
and the tests assert the individual and the file together. Under the mutation that takes the name
from the stem, three tests now fail.

**Two doc claims that were wrong about my own fixture**, both corrected against a run rather than
against the sentence:

- *"the difference is 40 bytes"* where the test asserts 30 — the fourth psp carries no census and
  costs nothing.
- *"the psp with no census is the fourth, so a pass that read every trailer but the last still
  fails"*. It does not: the byte count is 30 either way and only the row count notices. The
  census-less psp is now the **first**, which makes the two assertions catch different faults, and
  both numbers in the doc are measured: a pass that never reaches the last psp reads **20 bytes
  against the 30 asserted**, and a pass that judges every psp and then drops one leaves the byte
  count at 30 and is caught by the row count, **3 against 4**.

**And the reason given for carrying the path was a cohort the opener refuses.** Both doc comments
said *two directories of the same accession are one sample name and two files* — which
`refuse_a_sample_named_twice` rejects, about 460 lines below the second copy of that sentence. The
conclusion holds and the reason is now the true one: a sample name does not say which file it came
out of, and the path is what the user copies into `regenerate-census`.

**Smaller, all fixed:** `#[must_use]` on the pass, which returns the whole judgement and refuses
nothing if dropped; a `debug_assert_eq!` on the two vectors the accessor zips, since `zip` truncates
silently and what it would drop is one sample's row; *"one seek and ten bytes a sample"* softened to
*at most*, since a psp with no census costs neither; `fn judged` renamed to a clause, as it shadowed
three test bindings; the version-word fixture reused by the older test that still inlined it; the
cohort fixture taking `&[u8]` rather than `Vec<u8>`; and a duplicated import.

**Added on the range rule** (`CLAUDE.md`): a cohort of one. One stale psp is one row — not zero, and
not a refusal — which is the case a pass written as *report the samples after the first* would get
wrong.

### Ten mutations, all run, all caught

Each was applied to non-test code from a backup of the file, the module's tests run, and the file
restored and compared against the backup before the next.

| mutation | outcome |
|---|---|
| stop at the first psp whose census is not fresh | 4 tests fail |
| a psp that will not read reported as carrying no census | 2 fail |
| a psp that will not read dropped from the report instead of named | 1 fails |
| the rows sorted by sample rather than left in the order given | 2 fail |
| every row given the cohort's first path | 4 fail |
| each psp paired with another psp's path and name | 3 fail |
| the individual read off the file's name instead of out of the header | 3 fail — and **no assertion at all before the fixture was fixed**, every psp having been written to `<sample>.psp` |
| every psp judged and the last row then dropped | 5 fail |
| the last psp never read at all | the byte count fails, 20 against 30 |
| the whole trailer read where its front is enough | 3 fail, B2's own among them |

### What was measured

- **`cargo test --lib --bins --tests --all-features --no-fail-fast`: 6,705 lib tests pass** against
  the baseline's 6,682 — B1's four, B2's ten, B3's four and this step's five — with 20 of 21 targets
  green and the one pre-existing failure unchanged.
- **clippy: the same 11 errors of the same 5 kinds in the same 6 files as the baseline.**
- **`cargo check --all-targets --keep-going`: the same 4 examples.**
- **`cargo fmt --check`: the same 4 files.**

---

## Checkpoint B — what a run can now say about a cohort's censuses, and what it still cannot

**Milestone B is four steps and no behaviour a user can see yet.** What it builds is the judgement
the next two milestones spend: given a cohort of psps, *which samples need their census
regenerated, and why*. `estimate-parameters` starts asking at plan step C3; `regenerate-census`
at D2.

What holds after B4:

- **A psp's census is judged from the file's head**, for one seek and ten bytes — no census is
  decoded, so the answer is available before the run's reference is opened. It reaches two of the
  three causes spec §4.2 lists: no census at all, and a census of a format this build does not read.
  The third — a census built against a different set of loci — needs the reference read and the
  selection rebuilt, and stays where it is.
- **A cohort is judged whole**: one row a sample, in the order the run was given them, whatever any
  earlier psp turned out to be. A psp that will not read is a read failure in its own row, not a
  stale census.
- **A cohort whose psps were walked under different settings is refused when it is opened**, by
  whichever command opened it, naming both samples and the field that differs.

**Three things for the owner, none of them acted on.**

1. **The spec has a factual error in the sentence that argues for B3.** `psp_census_pair.md` §6 says
   `SegmentationInputs::first_difference` "is called only from its own tests". It is called from
   `PspVariantCaller::open`, and was before this branch — so a calling run already refused a cohort
   typed two ways, naming one sample and the run rather than the pair. The decision is unaffected:
   the check still belongs in the shared opener, which is what gives the fit and `regenerate-census`
   a refusal they had no way to make. The spec was left as it is.
2. **A census that decodes badly past its version word has no verdict.** Magic and version right,
   sections truncated: no cheap read tells it from a whole one, so it reaches the user as the census
   reader's own error rather than as *regenerate this*. Recorded in `CensusVerdict`'s own doc.
3. **`OpenPspCohort` lives in `psp_caller.rs`, which the module doc calls psp mode's *calling*
   stage**, and it now has four consumers outside that module — a fifth arrives at plan step C2,
   when `estimate-parameters` starts opening its cohort this way. The reviewer's suggestion is to
   lift the opened cohort into a module of its own at C2. It is a move, not a design change, and it
   is the owner's to time.
