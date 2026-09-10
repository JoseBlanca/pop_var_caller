# The census lives inside the psp — Milestone C: `estimate-parameters` over psps

**Date:** 2026-09-10
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone C
**Spec:** [psp_census_pair.md](../../ng/spec/psp_census_pair.md) §4, §4.1, §4.2, §4.3, §5, §6
**Branch:** `census-vs-psp-perf`, on top of `d3cfff1a`
**Milestone B's report:** [ng_psp_census_pair_milestone_b_2026-09-09.md](ng_psp_census_pair_milestone_b_2026-09-09.md)

---

## What this milestone is for

**Today a person running the parameters fit types the walk's settings a second time.**
`estimate-parameters` takes `--census` paths and five flags — `--min-copies`, `--min-period`,
`--max-period`, `--max-str-len`, `--min-purity` — that must equal what the psps were walked under,
and nothing checks that they do. A mistyped one does not fail: it rebuilds the run's selection of
loci over a different partition of the ground from the one the censuses hold, and what the person
sees, twenty seconds later and only sometimes, is the fit refusing the cohort as *built under
another selection*.

Milestone C takes those settings out of the user's hands: they come from the psp headers, which
recorded them at the walk. The command takes psps instead of census files, and a cohort whose
censuses are stale is refused in one report that names every stale sample rather than the first
one.

---

## The tree is not green, and every step is judged against this and not against green

Unchanged from Milestone B's baseline, and re-measured at `d3cfff1a` before C1 was written.

| gate | on the tree at `d3cfff1a` |
|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | **6,705 lib tests pass**, 15 ignored; 20 of 21 targets green; the one failure is `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` (`tests/ng_calling_loop_calls_genotypes.rs:1241`) |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | **11 errors of 5 kinds in 6 files** |
| `cargo check --all-targets --keep-going` | **4 examples** do not compile |
| `cargo fmt --check` | **4 files** |

The 6,705 is Milestone A's 6,682 plus Milestone B's 23 tests. **Every gate is compared as a set,
not as a count** — clippy's error kinds and their locations, the targets that will not compile, the
files `fmt` would rewrite — for the reason Milestone B's report gives: a check written from the
baseline's own failures cannot see a failure it has never seen before.

---

## C1 — the criteria from the header

**Committed:** see `git log` for `feat(ng): C1`.

### What it does

`estimate-parameters` used to build the run's segmentation from its own `--min-copies`,
`--min-period`, `--max-period`, `--max-str-len` and `--min-purity`. It now builds it from the
criteria the cohort's psps carry in their headers, which is what their walk recorded. **The five
flags are still on the command line and nothing reads them** — the plan asks for that split
deliberately, so that this step changes where one number comes from and nothing else; step C2
deletes them along with `--census`.

`run_ground::segments_over` is split in two at the line the spec's trap paragraph names. The flag
half converts the five flags into a `StrRepeatCriteria` and hands it on; the shared half —
`segments_cut_with` — refuses a catalog that is not there, checks the one that is against the
reference, cuts the ground and records what it was cut with. A psp-driven command enters at the
second, because a psp header carries the finished criteria and there is no conversion back to
flags.

### What the flags' help now says, and one refusal that changed

The five flags say **"Not read"**, and `--min-copies` says where the value comes from instead. A
flag whose help still read *"Give the values the psps were walked under"* while the code ignored it
is the state in which a stale sentence is believed.

One refusal is gone with them: a period range typed backwards (`--min-period 6 --max-period 2`) was
refused before anything was read and is now accepted in silence, because nothing converts those
flags any more. It returns at C2 as *no such flag*.

### The refusal that would have named a knob that does nothing

A catalog is built below every calling floor so a reader filters rather than re-scans; asking it for
tracts *below* what it holds is refused, because those rows were never written. That refusal names
the flag to move — and on the header-driven path there is no such flag, so a person pointing
`--catalog` at a catalog coarser than their walk would have been told, measured:

> `--min-copies asks for repeats the catalog /…/twenty-copies.repeats.parquet does not hold; raise
> it, or rebuild the catalog at lower floors`

— advice they cannot act on, since `--min-copies` is not read. `segments_cut_with` now takes
`CriteriaSource`, and on `ThePspHeaders` the refusal names the catalog instead, which is the thing
the reader can change. A test pins both halves: the message must name the catalog file and must not
name the flag.

### The oracle, and what it can and cannot show

**The fitted parameters file on the walked fixture cohort is byte-identical before and after: 20,741
bytes, same SHA-256.** That is the right check for *nothing else moved*, and it is all it is: the
fixture walks and fits with the same default criteria, so the flags and the headers carry the same
value there and no oracle of that shape can show which of the two the code read. What shows that is
the mutation table below — the first two entries.

### What the review found, and what was done

- **A refusal that names a flag the header-driven path does not have** — above. Fixed with
  `CriteriaSource`, and pinned by a test that builds a catalog holding tracts of twenty copies and
  up and points the fit at it.
- **Four statements in new prose that the code contradicts.** *"The ground and what it is cut with
  are the psps' own"* — they are the **first** psp's, and what makes one psp enough is a check
  elsewhere (`CohortCensusEvidence::new` refuses a cohort whose censuses disagree on their
  recording terms; at C2 the cohort opener does it directly). *"The conversion below checks the
  same thing again"* — the conversion checks nothing about the catalog; the other half of the split
  does. A doc listing *"finding the catalog"* as this function's work four lines above saying its
  path arrives resolved. And a test doc saying that before this step the second fit *cut the ground
  differently*, where in fact it would have been **refused** as built under another selection.
- **A claim more confident than the code supports.** The comment leaned on *the fit refuses a
  selection built under other criteria*. The fit compares a digest of the kept **generic**
  positions, so it catches a difference only where one of those moves; at tomato's 1-in-400 keep
  rate a criterion that retypes one short tract can pass it. That is an argument for taking the
  criteria from the psps, not a net to trust, and the comment now says so.
- **The catalog-path default was written out four times** — `<reference>.repeats.parquet` when
  `--catalog` is absent — with a fifth arriving at `regenerate-census`. One `catalog_path_for`, and
  `GroundRequest::catalog_path` calls it.
- **Two tests were satisfied by any successful fit.** The `--catalog` test now compares the file
  fitted from the moved catalog against the same cohort fitted with its catalog left where a run
  looks by default: one catalog read from two paths is one answer.
- **Nothing pinned that the ground comes from the cohort**, and no fixture here could: every walk
  covers whole contigs, where *the psps' ground* and *every base of the reference* are the same
  positions. The fixture takes a BED now, and a cohort walked over the first 400 bases of a
  600-base contig is fitted with no ground given.
- **Nothing pinned the refusal order** the new comment claims — a missing catalog before a
  backwards period range. The test named for that order calls the conversion directly and never
  meets a catalog. A test in `call_from_alignments` now sets both.
- Smaller: `let open = …` for an opened catalog renamed (the parameter took the name), the
  `# Errors` list reordered to the order the errors are raised, and the module doc now says the psp
  header is read for its criteria as well as its ground.

### Recorded, not fixed

- **Nothing compares the catalog the fit reads against the catalog the psps name.** Spec §6 says
  `--catalog` is "checked against the header's catalog digest as today"; there is no such check
  today. The header is in the cohort's own `SegmentationInputs` and
  `SegmentationInputs::first_difference` already names the field. **It belongs in C2's cohort
  opener** and is the last way a person can hand this command a wrong value and get numbers back.
- **The `Segmentation` this command builds is never read.** Its only use is
  `segmentation.inputs().repeat_tract_criteria`, which is the value handed in. The call is worth
  keeping for its refusals — a catalog that cannot serve the criteria is caught there — but it
  materialises every typed region of the analysed ground to read back a struct the caller already
  held, and that is the one cost of a run that grows with the genome rather than with the cohort.
  C2 rewrites this block and should either pass the criteria straight through or ask `run_ground`
  for the refusal without the segments.
- `generate-census` still takes its criteria from its own five flags, which is the shape
  `estimate-parameters` has just left. Milestone D replaces that command.

### Eight mutations, all run

Each applied to non-test code from a backup, the tests run, the file restored and compared against
the backup before the next.

**One of them had to be run twice, and the reason is worth keeping.** The first run of the
`--catalog` mutation reported *all tests pass* — because the review's fix had rewritten the line it
patches, so the edit matched nothing and the tree was never mutated. A mutation that does not apply
looks exactly like a mutation nothing catches. The harness asserts its match count; what hid it was
sending the apply step's output to `/dev/null` in the same command as the test run.

| mutation | outcome |
|---|---|
| the criteria taken from the five flags again, as before this step | 2 tests fail |
| the copy floor alone taken from `--min-copies` | 2 fail |
| `--catalog` ignored and the sibling path always read | 2 fail |
| the ground taken from the reference rather than from the cohort | 1 fails — the new ground test, and nothing else |
| the catalog that cannot serve the criteria refused naming `--min-copies` again | 1 fails |
| the catalog checked *after* the flags are converted | 1 fails — the new refusal-order test |
| the ground cut with `StrRepeatCriteria::default()` | 0 in this module, **4 elsewhere in the lib suite** |
| the criteria not recorded in `Segmentation::build` | 0 in this module, **5 elsewhere** |

The last two are the honest result and worth keeping: on this command the cut is genuinely unused —
only the record is read — so no test of `estimate-parameters` can see either, and what catches them
is `call-from-alignments` and `generate-psps`, whose runs do use the segments.

### What was measured

- **`cargo test --lib --bins --tests --all-features --no-fail-fast`: 6,711 lib tests pass** against
  the milestone's baseline of 6,705 — this step's six — with 20 of 21 targets green and the one
  pre-existing failure unchanged.
- **clippy: the same 11 errors of the same 5 kinds in the same 6 files.**
- **`cargo check --all-targets --keep-going`: the same 4 examples.**
- **`cargo fmt --check`: the same 4 files.**
- **The fitted file on the fixture cohort: 20,741 bytes, byte-identical to the pre-C1 run**, checked
  again on the finished tree after every review fix.
