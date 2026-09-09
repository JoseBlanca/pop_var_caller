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
