# ng STR observations — B2: the parameters file records what the run routed with

*2026-09-02. Step B2 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §2.3](../../ng/spec/run_ssr_observations.md) and adding
[`parameters_file.md`](../../ng/spec/parameters_file.md) §3.9. Branch `ng-ssr-observations`.*

## Plan

B1 made what counts as a repeat a knob. Two runs over the same reference and the same catalog
can now analyse different ground, and nothing a run wrote down said so. Documentation first —
a new §3.9 in the parameters-file spec settling the record's TOML spelling — then the record
itself.

## Assumptions

**The section is optional on read and always written.** The spec's ruling is *visible, never
blocking*, and a required section would block: `deny_unknown_fields` plus a required field
makes a file from an older build a parse error. `Option` is also §5's own idiom — absence is a
claim (*this file does not say*), never a sentinel for the defaults. Written into §3.9 as its
sixth row and pinned by a test.

**All eight axes of the criteria are recorded, not the five that have flags.** A record that
omitted the flank floor, the score floor and the bundling distance would not say what the run
actually asked the catalog for. The three are marked in the file's own comment as having no
flag, and the two that a difference in them implies — *rebuild the catalog* — is what the
comparison says by naming them in words rather than as flags.

**A file read in without a routing record and written out again gains one.** The written file
describes *this* run, and this run routed with something. The round-trip tests that assert
`of_run(file.project()) == file` therefore needed the fixture to carry a record; it now does.

**The comparison is reported through `eprintln!` beside the existing "replacing the parameters
file already at …" note.** That is the channel this command already uses for a fact an operator
must see that stops nothing, and it needed no new surfacing decision. Nothing in the spec named
a channel.

## Changes made

| file | change |
|---|---|
| `doc/devel/ng/spec/parameters_file.md` | §3.9, and a sixth row on §5's absent/zero/default table |
| `parameters_file/mod.rs` | `RepeatRouting`; `ParametersFile::repeat_routing: Option<RepeatRouting>`; the every-shape fixture carries one |
| `parameters_file/bindings.rs` | `RepeatRouting::of`; `ParametersFile::routing_disagreement`; the demotion answers for the new field |
| `parameters_file/from_run_parameters.rs` | `of_run` takes the run's criteria and records them |
| `parameters_file/to_toml.rs` | the section, its comment block, and `a_toml_purity` |
| `call_from_alignments.rs` | the criteria come from the segmentation; a supplied file that differs is a `note:` line |
| `testdata/*.toml` | both goldens regenerated through their own `--ignored` regenerators |

**`a_toml_purity` exists because widening the floor lied about it.** The purity floor is an
`f32`; formatting it as the `f64` it becomes spells the shipped 0.8 as `0.800000011920929` —
every digit true of the `f64` and none of them anything a person typed. Formatting the `f32`
gives `0.8`, and it survives the trip through TOML's `f64` because an `f64` carries more than
twice an `f32`'s precision, so narrowing back cannot land on a different `f32`.

**The record is read off the segmentation, not rebuilt from the flags.** Rebuilding would let
the file say one thing while the catalog had been asked another; `Segmentation` already keeps
the criteria the reader asked with, for the compatibility checks, so there is one value.

## Tests added

Six, and the last is the one that could fail silently.

- `the_record_carries_every_axis_the_catalog_was_asked_on` — all eight, each read back
  separately. A record dropping one would let two runs that routed differently write identical
  files.
- `the_same_routing_is_no_disagreement` — the control.
- `each_axis_that_differs_names_itself` — one case per axis, so a comparison that stopped at the
  first field, or compared the record against itself, fails.
- `a_file_with_no_routing_record_makes_no_claim_to_disagree_with` — §5's rule: a run must not be
  told a file disagrees when the file said nothing.
- `a_routing_difference_leaves_every_number_as_warranted_as_it_was` — the comparison is a
  question and not a step; the file is unchanged by being asked.
- `the_written_parameters_file_records_what_this_run_counted_as_a_repeat` — end to end through
  `run_call_from_alignments`, with **all five flags moved off every default**, so a writer that
  recorded the catalog's floors, or the calling defaults, or the catalog file's own header,
  fails.

**Mutation-tested**: making `of_run` record the calling defaults regardless of what the run was
given — the exact silent failure this step risks — fails the end-to-end test and none of the
other five. Restored from a backup and the mutation's absence checked by grep before committing.

## Validation

In the dev container:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --tests --examples --all-features --no-fail-fast` — **5,926 passed, 0
  failed, 14 ignored** in the library suite (5,920 before this step); every integration target
  green.
- Both goldens regenerated through their own `--ignored` regenerator and the diff read: one new
  section each, nothing else moved.
- Unchanged by this step and still red: the three locus-dump tests and the psp writer bench,
  recorded in `PROJECT_STATUS.md`.

## Tradeoffs and follow-ups

- **The fit's own routing criteria are still only a digest.** They sit inside
  `[fitted_from.census]` as `STR routing criteria`, which answers an equality and cannot be
  read by eye. So a run comparing itself against *its fit* still cannot name the axis; what it
  can do is compare against a parameters file written by another ng run, which is what
  `routing_disagreement` does. Making the census term legible is the census format's, not this
  step's.
- **`FORMAT_VERSION` is unchanged at 1.** The section is additive and optional, so a file from
  either build reads in the other; a version bump would imply a migration policy that spec §11
  defers until there is a file to migrate.
