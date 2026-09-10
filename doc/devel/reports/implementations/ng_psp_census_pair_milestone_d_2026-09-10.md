# The census lives inside the psp — Milestone D: `regenerate-census`, the repair

**Date:** 2026-09-10
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone D
**Spec:** [psp_census_pair.md](../../ng/spec/psp_census_pair.md) §3.1, §6, §8, §11
**Branch:** `census-vs-psp-perf`, on top of `93fb1ac3`
**Milestone C's report:** [ng_psp_census_pair_milestone_c_2026-09-10.md](ng_psp_census_pair_milestone_c_2026-09-10.md)

---

## What this milestone is for

**A fit now refuses a cohort and names the samples whose censuses it cannot use, and tells the
person to run a command.** Milestone C built that refusal; the command it names did not exist.
Today's `generate-census` writes a census *file* beside a psp, and nothing reads those any more —
so following the refusal would cost a quarter of an hour a sample and leave every psp exactly as
stale.

Milestone D turns that command into the repair: psps in, each one's trailer replaced with a fresh
census, nothing else in the file rewritten. D1 is the rename and the write; D2 skips the psps that
need nothing; D3 is a run stopped part-way; D4 is the whole-file oracle.

**The baseline every step is judged against is the one Milestone C's report records** — the tree is
red in four ways, none of them this plan's — re-measured at `93fb1ac3` in each step's gate table.

---

## D1 — the rename, and the write into the psp

**Committed:** see `git log` for `feat(ng): D1`.

### What it does

`generate-census` is `regenerate-census`. The module, its tests and its `SUBCOMMAND` moved with
`git mv`, so the history follows the file, and what a person can type is now the whole of spec §6's
list: `--psp`, `--reference`, `--catalog`.

Per sample: `census_from_psp` reads the psp's records and builds the census, `write_census` encodes
it **with no pileup identity** — a census that *is* its psp's trailer has nothing left to pair
wrongly with (spec §3.1) — and `replace_trailer` puts it where the walk had put one. The header,
the blocks and the index are the bytes they were.

**The cohort is opened to be agreed with, then closed before the first psp is rewritten.** Opening
it is what refuses files walked over different ground, under a different catalog or different
criteria, or sharing an `@RG ID` (step B3), and what says what the ground and the criteria are.
Holding a thousand psps open to rewrite them one at a time would spend the memory psp mode exists
to save — and each one is about to be truncated at its trailer.

### Three things absorbed beyond the plan's own sentence for this step

The plan's D1 line says `--output-dir` and `--force` go. Two more changes came with it, and the
step's review was asked to challenge all three.

1. **The five criteria flags went too.** Spec §6 is explicit — this command takes `--psp`,
   `--reference` and `--catalog` and nothing else — and the criteria come from the psp headers,
   which is what C1 and C2 did for `estimate-parameters`. Left in, they would be values a person
   retypes from a walk that already recorded them, and a mistyped one rebuilds a census under
   settings they think they did not change.
2. **The reference and the catalog are checked against the psps before anything is rewritten.**
   This command *writes* the census, so a wrong file here is worse than at the fit: it produces a
   census whose recorded settings name another reference or another catalog, which every later fit
   refuses — after this run has spent the rebuild. C5's review measured that hole in the command
   this step replaces: it can write a census recording one reference into a psp whose header names
   another. Both refusals name the file and say to run again with the psps' own; the reference
   refusal names the FASTA the psps record, out of their headers.
3. **The catalog comparison is the catalog header's own method now**
   (`RepeatCatalogHeader::first_difference`), so this command and `estimate-parameters` name the
   same difference in the same words. It has its own tests, and they pin the clause *order*: the
   reference first, because a catalog built on another assembly makes every other comparison
   meaningless.

### What the old command refused, and where each refusal is now

| refusal | now |
|---|---|
| the reference will not read | kept |
| the ground, the selection, the catalog | kept |
| a `--psp` directory that will not list, or holds no psp | kept, through the rule the three psp-taking commands share (`psps_named`) |
| the cohort's files do not agree | kept, and it is the same opener `call-from-psps` and the fit use |
| a sample name that cannot be a file name | **gone with the file it was about.** It existed because `<sample>.census` was a path built from `@RG SM`, which is free header text; nothing derives a path from a sample name now |
| a census already at the output path, and `--force` to replace it | **gone.** There is no output path: the census goes into the psp, and replacing it is the command |
| — | **added:** a reference the psps were not walked against; a catalog they were not walked with; a psp that will not take the write, which says whether that file was left as it was or cut back to its trailer |

### The parity oracle in its new shape

The old oracle compared a census *file* this command wrote with the census in the psp's trailer.
There is no file now, so it compares the trailer the walk wrote with the trailer this command writes
in its place — on the fixture with a repeat tract and three libraries, and on the one whose second
sample has no reads at all.

**What it still covers, and it is the part that matters:** the selection is derived twice. The walk
builds its segmentation and its `CensusPlan` from its own flags, this command builds them from the
psp headers, and a divergence about how a plan is *built* changes which loci are kept and therefore
the bytes — measured by mutation, seeding the selection differently fails both parity tests.

**What it did not cover, and the step's review caught both.** A comparison of a psp's trailer before
and after the run holds just as well when nothing is written at all: with the write removed, all
fourteen tests passed. And the record count the report prints was compared only against the field
the report was built from. Both are pinned now — one test empties a psp's trailer so the bytes it
asserts can only have come from this run, and the line's count is asserted against the psp's own.
The census *file's* layout and the pileup identity, which a reader might expect in this list, are
covered elsewhere: by `census_file.rs`'s round-trip tests and by a test of this command's own.

### A wrong number of my own, caught by its own test

I wrote that the plain fixture's two samples declare three libraries between them. They declare
two, one each; the assertion failed on its first run and now says the measured number. It is the
failure the plan-driven skill names as the most reliable defect this loop produces — a claim about
my own fixture, recalled rather than measured.

### What was measured

- **6,731 lib tests pass**, against 6,728 on the tree C5 committed: the command's thirteen tests
  became fourteen, and the catalog-header comparison brought two of its own.
- `cargo check --all-targets --keep-going` fails on exactly the baseline's four examples, so
  **nothing that compiles into the gate names the command that was renamed**. Four things outside it
  still do: two scripts and one example, which are plan step E1's own work
  (`scripts/ng_fit_stage_end_to_end.sh` is now broken rather than out of date — it passes
  `--output-dir` to a subcommand that no longer exists — plus `scripts/ng_census_route_cost.sh` and
  `examples/ng_census_route_cost.rs`), and three doc comments in `src/`, which the step's review
  found and this step's follow-up commit fixes.
- Full gates are in the commit message, compared as sets against the milestone baseline.
