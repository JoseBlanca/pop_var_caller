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
