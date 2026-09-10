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

---

## C2 — `--psp` in; `--census` and the five criteria flags out

**Committed:** see `git log` for `feat(ng): C2`.

### What it does

`estimate-parameters` takes psps. `--census` is gone, and so are `--min-copies`, `--min-period`,
`--max-period`, `--max-str-len` and `--min-purity`: there is nothing left on this command line that
a person could type differently from the walk that wrote the files. What a run says now is the
reference, the catalog, the psps and where to write.

Three pieces:

- **The cohort is opened by the opener every psp-taking command shares.** `OpenPspCohort::open`
  reads each header and refuses a set that was not walked as one — two files naming one sample,
  different ground, a different catalog, different repeat-tract criteria (step B3) — and no block
  is decoded by any of it.
- **The censuses come out of the psps' trailers**, through A3's reader, which takes a path and an
  extent. Each census's header and directory are decoded at open; a section is read when the fit
  asks for it.
- **`open_census_cohort` and `CensusInCohort` are gone**, and with them the class of refusal that
  existed because a census was a separate file.

### What the old opener refused, and where each refusal lives now

| the old refusal | now |
|---|---|
| nothing named | `RunError::NoPsps` in the psp opener, and clap's own `required` |
| a file that is not a census | the same variant, keyed on the psp that carries it |
| a psp that will not read | `PspReader::open`, which validates the footer, the index bounds and the header, where the old check read the header alone |
| two samples walked over different ground | `the_settings_every_file_agrees_on`, which compares the catalog and the criteria as well, and names both samples and the field |
| a census with no psp beside it | gone: a trailer cannot be separated from its file |
| a census that names no psp | gone, same reason |
| a census built from another psp | gone, same reason |

**One input is newly unrefused, and it cannot arise today**: a psp whose trailer holds *another
sample's* census. Only the walk writes a trailer, and it writes each sample's own. `regenerate-census`
(plan step D1) will be the second writer, and that is where the check belongs — the review raised it
and it is recorded below rather than added here, because a guard with no reachable failure is a
guard nothing can test.

### What it costs, which is not nothing

**The fit now holds every psp open for the whole run.** Before, it read one header per sample and
dropped it. An open psp keeps its block index — 24 bytes an entry and about 14,000 entries for a
whole genome, so **about 336 kB a sample** by the format's own arithmetic
(`psp/reader.rs:82-86`), which is roughly 340 MB at a thousand samples. That is the shape spec §5
asks for and the same shape `call-from-psps` already has; it is named here because it is the fit's
first per-sample resident cost that grows with the *genome* rather than with the census, and
because a fit-side budget is explicitly out of this plan's scope.

Bytes read a sample are otherwise as before: the census reader's head read, at most a megabyte,
out of which a few hundred bytes of header and directory are decoded.

### The oracle

**The parameters file fitted over the fixture cohort is byte-identical to the one fitted before
C1: 20,741 bytes.** It is now fitted from censuses read out of psp trailers where it was fitted
from census files, so this is A1–A4's byte-for-byte parity carried through the whole fit — the
walk's census, the census `generate-census` wrote, and the file the fit produces from either.

### What the review found, and what was done

**Two blockers, and the second is the one worth reading.**

1. **A rustdoc link to the function this step deleted**, and broken intra-doc links are denied, so
   `cargo doc` would have failed on it. (While fixing it: `cargo doc` **is already red on this
   tree** — 40 unresolved links, none of them from this branch. It is a fifth gate nobody has been
   running, and it is not this plan's to fix. Raised at Checkpoint C.)
2. **The new laziness test could not fail.** It asserted that reading a cohort's censuses costs
   fewer bytes than the trailers carry, using the census reader's own counter — which counts
   *section* reads and not the head read that opening does. So the number was zero under every
   implementation, including one that took each trailer whole, and the assertion reduced to *the
   psps have a non-empty trailer*. **The instrument that does discriminate is the psp's own**
   (`trailer_bytes_read`), and the replacement asserts three things: no trailer taken whole, no
   section decoded, and then — asking for a section afterwards — that the read reached the disk,
   which is the only way to tell a lazy census from one already in memory. The mutation that reads
   each trailer whole and decodes it resident now fails it.

**Then eight smaller things, four of them prose the code contradicted**: "six tests went" where
five did and only three were about the pairing; "a few hundred bytes" for a head read of up to a
megabyte; "the five flags are still on the command line, and C2 removes them" in the step that
removes them; and a stale error doc naming the deleted opener. Also: `&mut OpenPspCohort` where the
function seeks nothing, now a shared borrow with its own accessor; two variables called `censuses`
and `open` holding psps and evidence; `CohortRefusal`'s three causes described as two; and a
refusal test that asserted only the outer error variant, now naming the psp.

**And the `--psp` listing rule, written out twice.** `call-from-psps` and `estimate-parameters` had
the same thirty lines — expand a directory, keep the psps directly inside it, sort by name, refuse
an empty one — differing only in which error type they built, with a third copy due at
`regenerate-census`. It is one module now (`psp_inputs`), with five tests of its own, and each
command dresses the refusal in its own words. The rule is not obvious enough to be safely written
three times: the sort **is** the run's sample order, and it reaches the VCF's columns.

### Five mutations, all run, all caught

| mutation | outcome |
|---|---|
| each trailer taken whole and decoded resident | 1 test fails — the rewritten laziness test |
| a census that will not read skipped instead of refused | 2 fail |
| each psp paired with another psp's path | 6 fail |
| the evidence's sample order reversed | 1 fails |
| the read-group table's file column taking the wrong path | 3 fail |

The last was the review's prediction of a survivor — nothing read that column — and it is now
pinned by one assertion, because the column exists so that a message can name a file and a message
naming the wrong file is worse than one naming none.

### What was measured

- **6,712 lib tests pass** against the milestone baseline's 6,705: C1's six, less the four tests
  whose subject this step deleted, plus the shared listing module's five.
- **clippy, `cargo check`, `cargo fmt --check`: the same sets as the baseline.**
- **The fitted file: 20,741 bytes, byte-identical to the pre-C1 run**, re-checked after the review
  fixes.

### Recorded, not fixed

- **A trailer holding another sample's census is not refused.** Free to check — both names are in
  hand — but unreachable until `regenerate-census` exists, so it belongs with that command (D1).
- **`cargo doc` is red on this tree**: 40 unresolved intra-doc links, none from this branch.
- The end-to-end script's step 3 now passes `--psp`; the rest of that script — the
  `generate-census` step it no longer feeds — is plan step E1's.

---

## C3 — the refusal, before the reference

**Committed:** see `git log` for `feat(ng): C3`.

### What it does

A run whose psps do not all carry a census this build reads now stops with a message that names
every one of them, and ends with the command that rebuilds them:

```
2 of this cohort's 4 samples cannot be fitted as they stand:
  bravo (/psps/bravo.psp) carries no census
  delta (/psps/delta.psp) carries a census built by an older version of this program: it is
    version 1 and this build reads version 2
2 other samples carry a census this build reads.
Rebuild the 2 samples whose census is named above, then run this fit again:
  regenerate-census --reference ref.fa --catalog ref.fa.repeats.parquet --psp /psps
```

**Every psp is judged first, and the run stops after all of them, not at the first.** Rebuilding
one census is a quarter of an hour a sample (spec §2, §4), so a refusal that named them one at a
time would cost that wait once a stale sample, in series, to learn a job that fits in one message.

**And the judgement happens before the reference is opened**, which is what makes the refusal
immediate: judging a psp is one seek and ten bytes a sample, where reading a human reference is
minutes. A test pins it by pointing `--reference` at a path that does not exist and asserting the
run still answers with the report.

**The stale lines are grouped by cause** (spec §4.1), keeping the run's sample order inside each
group, so *these three were written by an older build* reads as one job.

**A psp that cannot be read is reported apart, and the command is not offered for it.** Rebuilding
the census of a truncated file would not mend it. The type carries the command as an `Option`, set
only when something is actually stale, so a report cannot tell someone to rebuild what rebuilding
will not fix.

### The command it names does not exist yet

The report says `regenerate-census`, which arrives at plan step D1. **Naming today's
`generate-census` instead would be worse**: that writes a census *file* beside the psp, which this
fit no longer reads, so a person who followed it would spend the wait and find their psp exactly as
stale. The name is a constant with one place to point at D1's real subcommand, and a test asserts
it is not the old command's name. At D1 that test becomes what it should be — a parse of the name
against clap, which cannot be written today because it would fail.

### What the review found, and what was done

No blockers. **The message was the review's subject, judged as what a person reads at 2am when a
60-sample fit stops**, and three of its findings were about that:

- **It inflected one of its three number-bearing phrases.** At one stale sample it said *cannot be
  fitted as they stand* and *Rebuild them*; at a cohort of one — the low end of the range this
  caller is built for — *1 of this cohort's 1 samples*. Every phrase inflects now, and a test reads
  the singular case whole.
- **"Rebuild them" included the psps that cannot be rebuilt.** In a cohort with both faults, the
  unreadable rows printed directly above that line. It names the set now — *rebuild the 3 samples
  whose census is named above* — and says in its own line that the unreadable ones will not be
  mended by rebuilding. A test covers the mixed case, which had none.
- **The unreadable row printed its path twice and dropped the real cause.** It stored the
  outermost message where the fault — *unexpected end of file* — hangs off the source. It walks the
  chain now, and the fixture's error carries the row's own path so the duplication would show.

**The only arithmetic in the message was asserted by nothing**: the first line, which is what a
reader acts on before anything else. Changing *of this cohort's N samples* to something else broke
no test. It is asserted whole now, and the mutation fails three.

**Smaller:** three-field tuples where two row structs belong; an error variant named
`CensusesNeedRegenerating` for a report that also covers files nothing can regenerate, now
`CohortCannotBeFitted`; the command line built on every run including the successful one, now a
closure called only on the refusal path; a printed path with a space in it, now quoted, because the
whole value of that line is that it is pasted; and two doc comments that still described a census
file beside its psp.

**And one claim of mine that was false**: *every psp is judged before anything else is read*. By
that point the cohort opener has read every header, every footer and every block index — a few
hundred megabytes at a thousand samples. What is true, and load-bearing, is that it happens before
the **reference** is opened.

### Eleven mutations, all run

| mutation | outcome |
|---|---|
| the run stops at the first stale psp | 2 report tests and 2 command tests fail |
| the reference read before the judgement | 1 fails — the ordering test |
| the header's arithmetic changed | 3 fail |
| the sample-name column dropped | 3 and 3 fail |
| the cause dropped from each line | 2 and 3 fail |
| a psp that will not read reported as one to regenerate | 2 fail |
| the stale lines left ungrouped | 2 fail |
| the cause chain not walked for an unreadable psp | 2 fail |
| `--psp` dropped from the printed command | 2 fail |
| `--catalog` dropped from the printed command | 2 fail |
| the old command name printed | 1 fails — the test added for exactly it |

The last was the review's predicted survivor: the other assertions compare against the constant, so
changing the constant changes both sides. What can be pinned today is that the name is **not**
`generate-census`, and that is what the new test says.

### Recorded, not fixed

- **A psp that will not *open* never reaches the report.** `OpenPspCohort::open` refuses the cohort
  at the first file it cannot read, naming that one, so spec §4.1's *every sample is examined* does
  not hold for that fault — only for censuses. Making the opener collect its refusals the way B4
  collects verdicts is a change to the opener every psp-taking command shares, so it is the owner's
  call rather than this step's. The report's own doc says which fault it covers.
- **`CensusVerdict`'s doc promised grouping by cause and the first draft did not group.** It groups
  now, so the promise holds — but spec §4.1's wording and that paragraph should be read together at
  the checkpoint, because *grouped by cause* and *in the run's sample order* are two different
  reports and the code now does both, one inside the other.
- **Plan step D1's task list does not mention the constant** this step introduces, and `grep
  generate-census` will not find it, because it spells the new name.

---

## C4 — the fit's own refusal says what to run

**Committed:** see `git log` for `feat(ng): C4`.

### What it does

A census can be stale for three reasons (spec §4.2). Two are visible from the psp's head and are
C3's report. **The third — a census written against a different set of positions from the one this
run rebuilds — cannot be seen without reading the reference and rebuilding the selection**, so it
arrives at the end of that work, from the fit, as `CohortFitError::AnotherSelection`. Its message
explained what had happened and stopped. It now ends with what to do:

> …if the psps were walked against another reference or catalog, fit with those, which regenerates
> nothing; if not, this build chooses census positions differently from the one that wrote them, so
> regenerate them with regenerate-census and fit again

**Two faults end in this refusal, with opposite fixes, and the free one comes first.** A run
pointed at the wrong reference or catalog needs the right files; regenerating would rebuild every
census against the wrong reference — a quarter of an hour a sample at whole-genome scale — and the
next fit would be refused again. A build that chooses positions differently needs its censuses
regenerated, and no reference helps. The first draft offered both and led with the costly one; the
review caught it, and a test now asserts the order.

**The command's name is one constant**, beside the freshness judgement, read by both places that
tell a user to regenerate: this refusal and C3's report.

### What the command-level test found

**No test reached this refusal through the command**, so a command that turned it into any other
error passed everything. The test that now does rebuilds each psp's census under a selection with a
different position budget and runs `estimate-parameters`.

**Its first version asked for half the shipped budget, and the fit accepted the cohort without a
word.** The fixture's contig holds fewer ordinary positions than either budget, so both selections
kept the same set; the fit compares the kept set and not the terms it was chosen under, and every
census's recorded terms said the budget differed. The test uses a budget of three, which changes the
set. **The silent acceptance is Checkpoint C's first item.**

### Mutations

| mutation | outcome |
|---|---|
| the constant set to the old command's name | 2 tests fail |
| the fit's refusal swallowed into another error at the command | 2 fail |
| the instruction negated | 4 fail |
| the costly fix offered first | 2 fail |

### What was measured

Gates on the final tree are in the commit message: the lib suite gains the one command-level test,
and clippy, `cargo check` and `cargo fmt --check` are compared as sets against the baseline.

---

## Checkpoint C — step 2 takes psps, and what it can still be told wrongly

**Milestone C is four steps, and after it a user runs the parameters fit by naming psps.** The
command takes `--reference`, `--catalog`, `--psp`, `--output`, `--force`, `--ploidy` and
`--inbreeding`; nothing about what a repeat is, and no census file. What holds:

- **The repeat criteria and the ground come out of the psp headers**, and a cohort whose files
  disagree about either is refused when it is opened, naming the two samples and the field.
- **A cohort with stale censuses is refused before the reference is read**, in one message that
  names every stale sample, its file and its cause, and ends with the command that rebuilds them.
- **The fit's own refusal for a census built against another selection says what to run.**
- **The parameters file fitted over the fixture cohort is byte-identical to the one fitted before
  this milestone — 20,741 bytes** — now read out of psp trailers where it was read from census
  files.

### Four things for the owner

1. **`--reference` and `--catalog` are still values a person can type wrongly and get numbers
   back.** Each census records a digest of the reference it was built against and of the catalog's
   build settings, among seven selection terms. Those are compared sample against sample, and never
   against the selection this run rebuilds. What the fit does compare is a digest of the kept
   *generic* positions, which catches a difference only where it moves one of them — and on tomato
   about 1 position in 400 is kept (spec §2), so a catalog that retypes a few short tracts can pass
   it. **Spec §6 says `--catalog` is "checked against the header's catalog digest as today"; no such
   check exists, today or before this branch.** And it is not hypothetical: writing C4's test, a
   cohort whose censuses were rebuilt under **half the shipped position budget** fitted without a
   word on the fixture, because its contig holds fewer ordinary positions than either budget, so
   both selections kept the same set — while each census's recorded terms said the budget
   differed. Nor does `--reference` get compared with the psp headers, which `call-from-psps`
   already does with an existing function. The fix is small and the pieces are in hand: the
   rebuilt plan carries its selection terms, each census carries their digest, and a comparison that
   names the first differing field already exists. **Recommendation: add it as one more step before
   Milestone D**, refusing with the field's name before anything is fitted, with the positions
   digest kept as a backstop.
2. **One more stale-census refusal names no action.** A cohort walked by two builds with different
   selection constants passes the freshness judgement — every census has this build's format — and
   is then refused as *samples A and B disagree on selection seed; they did not record the same
   thing*: two samples named, nothing to do. Spec §1's goal is that a stale census stops step 2
   with every such sample, the reason and the command. **Recommendation: fold it into the same
   step as item 1** — it is the same comparison of recorded selection terms, run across samples
   instead of against the run, and its answer belongs in C3's report.
3. **A psp that will not open stops the cohort at the first one.** Spec §4.1's *every sample is
   examined before anything is refused* holds for stale censuses — C3's report — and not for a file
   that is truncated or missing, which the shared cohort opener refuses one at a time.
   **Recommendation: leave it.** An unopenable psp is damage rather than a regeneration job, the
   refusal names the file, and changing it means changing the opener three commands share.
4. **`cargo doc` is red on this tree, and it is not in any baseline.** Forty unresolved intra-doc
   links, none from this branch; C2's review found it because broken links are denied and C2 had
   just added a forty-first. **Recommendation: a separate clean-up, outside this plan**, and adding
   `cargo doc` to the gate set from then on.

### Two things to know

- **The fit now holds every psp open for the whole run.** An open psp keeps its block index — about
  336 kB a sample at whole-genome scale by the format's own arithmetic, roughly 340 MB at a thousand
  samples — where the fit used to read one header a sample and drop it. Spec §5 asks for that shape.
- **The command the refusals name does not exist until Milestone D.** Both messages say
  `regenerate-census`, from one constant. Plan step D1's task list does not mention it; D1 has to
  point it at the new subcommand's own name and replace the test that pins it with a parse against
  the command line.

---

## C5 — a census recorded under other settings is refused, naming the setting

**Committed:** see `git log` for `feat(ng): C5`. **Added to the plan at Checkpoint C by the owner's
ruling of 2026-09-10** — *"that's what we should do, stop and report the problem"* — and it closes
the first two of Checkpoint C's four items.

### What it does

`estimate-parameters` now makes five checks before it fits, and the order is the design:

1. **the psps' heads**, before the reference is read — a missing census, an older format, a trailer
   that is not a census (C3, unchanged);
2. **each census read out of its psp**, and **not yet assembled into a cohort**;
3. **`--reference` against every psp header** — the comparison `call-from-psps` has always made
   (`refuse_a_file_against_another_reference`), which this command never made;
4. **the catalog's header against the header the psps recorded**, before a segment is cut — spec
   §6 claimed this check existed and it did not;
5. **the selection rebuilt, and each census's twelve recorded settings against the run's own** —
   one verdict a sample, every stale one named, as rows of C3's report.

Only then are the censuses assembled into one cohort, and the fit's comparison of the kept
positions is left behind all of it as a backstop.

### The two faults have opposite fixes, and the order is what lets each message be unconditional

C4's refusal had to offer both ways out — *fit with the other files* and *regenerate* — because it
could not tell which fault it had met. Each is wrong for the other fault: a census regenerated
against the wrong reference costs a quarter of an hour a sample at whole-genome scale (spec §8) and
is refused again by the next fit with the right one.

Checks 3 and 4 catch the wrong-file fault, name the file, and say to run again with the psps' own —
which rebuilds nothing. What survives them is a census that does not match its own psp, or one
written by a build that chooses positions differently, and regenerating is the fix for both. So
each refusal now names one fix and does not hedge.

### The deviation from the plan's text, and the chain that justifies it

The plan asks check 5 to name the fix each setting calls for: a different reference or catalog
means *rerun with the files the psps were walked against*, anything else means *regenerate*. **As
built, check 5 always says regenerate**, because by the time it runs the run's reference and
catalog are provably the psps' own. The chain, each link read in the code and confirmed by the
step's review:

1. this command reads its reference from a FASTA, so the run's whole-assembly digest is always
   present; a `.fai`-only read cannot reach check 5 at all, because `ReferenceDigest::of` fails
   first;
2. a psp that carries a census was walked by a run that got past that same call, so its header's
   digest is present too;
3. with both present, check 3 compares them exactly — so after it, the run's reference has the
   psps' bases;
4. `open_catalog` proves the run's catalog carries the run's reference digest, and check 4 proves
   its whole header equals the psps'. `CatalogBuildSettings` — what the census records about the
   catalog — draws on three of those header fields, all compared;
5. the remaining four selection values come from the psps (the analysed ground and the routing
   criteria, both forced equal across the cohort when it was opened) or are this build's constants
   (the seed, the position budget, the per-stratum cap).

**And the case the plan worried about does occur, with the other answer.** Today's
`generate-census` makes no reference check, so it can write a census recording reference X into a
psp whose header says Y. Pointed at Y, check 5 names *reference digest* and says regenerate — which
is right, because regenerating rebuilds that census from the psp against Y. Pointed at X, check 3
refuses first.

### What a person reads

A cohort whose censuses were rebuilt under half the shipped position budget, which fitted without a
word at C4:

```
2 of this cohort's 2 samples cannot be fitted as they stand:
  one (…/psps/one.psp) carries a census recorded under settings this run does not use (the first that differs: generic target position count)
  two (…/psps/two.psp) carries a census recorded under settings this run does not use (the first that differs: generic target position count)
Rebuild the 2 samples whose census is named above, then run this fit again:
  regenerate-census --reference …/ref.fa --catalog …/ref.fa.repeats.parquet --psp …/psps
```

A wrong `--reference`, which used to be reported as a catalog built on another reference — blaming
the one file that was right:

```
…/another-build.fa is not the reference these psps were walked against; they name ref.fa, so run
this fit again with that one, which regenerates nothing: the psp for sample one was written
against a different reference from this run's: it was walked against the assembly whose checksum
is … and this run's reference is …
```

The psps' own name for the reference comes out of their headers, which record the FASTA's basename
for exactly this purpose. The catalog refusal cannot do the same, and says why in its doc: a
catalog's header holds no path and no name, so that message names the file this run read and what
about it is not theirs.

### What C5 makes unreachable, and what it delays

**Unreachable, and that was the point.** A cohort walked by two builds used to be refused as
*samples one and two disagree on selection seed; they did not record the same thing* — two samples
named and nothing to do. It cannot be reached from this command any more: if check 5 finds nothing,
every census records the run's settings, so no two of them can disagree. The plan asked for that
refusal to say what to do, preferably as rows of C3's report; it is now those rows.

**Delayed.** Assembling the cohort moved from before the reference read to after it, so the three
refusals that are about damage inside a census — two samples claiming one read group, two declaring
one `@RG ID`, a section for a read group the census does not declare — now arrive after the
reference, the catalog and the selection rebuild, which spec §4.2 measures at 4 to 19 s on the
tomato fixture and is minutes on a human reference. None of them is a mistyped command: the
ordinary shared-`@RG ID` case is still refused from the headers when the cohort is opened, before
anything is read.

### The run's own copy of what a census records, and how it is held in place

Check 5 needs the run to say what a census written under its own plan would record.
`CensusPlan::recording_terms` builds exactly that — the selection digest, the kept-position digest
in the writer's own order, the per-stratum counts, the two caps and the census's depth ladder —
from the same pieces `writer_for` hands the writer.

**It is not merely inspected, it is pinned by the ordinary case**: any divergence in any of the
twelve makes *every* fresh cohort stale. Two mutations measured that — a different depth ladder,
and the kept positions digested in reverse — and each failed 7 tests, including every test that
fits the walked fixture.

### What the review found, and what was done

One read-only agent over the grouped categories, forbidden to edit or build. **No blockers.** It
confirmed both structural claims above — the fix-split chain and `recording_terms` value by value —
and found five should-fix items, all fixed here:

- **the reference refusal told the user to switch references without saying to which.** The psp
  header records the FASTA's basename for this; it is now in the message, and asserted;
- **the catalog comparison's first clause could not fire** — both catalogs are checked against the
  run's reference as they open — **while the half of it that could said something false**: a
  catalog built from the same bases wrapped at another line width would have been reported as built
  on another reference. The dead half is now a `debug_assert` recording why it cannot differ, and
  the live half says *its contig table is not theirs*;
- **`CohortFitError::AnotherSelection`'s doc contradicted itself** — one paragraph said the cause
  arrives at the fit *rather than* in the report, eighteen lines above the new paragraph saying the
  command reports it. Rewritten, in the variant's doc and in its test's comment;
- **a comment claimed such a cohort is "refused twice, once for each"**, which spec §8 forbids:
  `regenerate-census` skips a psp only when all three causes are ruled out, so the one command line
  the first refusal prints rebuilds both psps and the next fit succeeds. Corrected;
- **the documented panic was stricter than the assertion**: the doc said a census paired with the
  wrong psp panics, and the code checked only the count. The pairing is now asserted by sample
  name, with a test that panics on a reversed list — and it earned its place immediately, by
  failing its own sibling test, whose fixture censuses all carry the fixture's sample name where
  the psps are named delta, alpha and charlie.

Ten minor items were taken too: four sentences of prose the diff had left describing the code it
replaced, the `# Errors` lists on the two halves of the `run_ground` split, a note that the
catalog's flag-naming path is now defensive for this command, and a note that the composed
census reader has no non-test caller left.

**Not taken, and recorded for the owner:** spec §4.2's third row and the plan's C5 text describe
what the code now treats as the backstop, and the plan says *seven* settings where the code
compares twelve. Widening to twelve is what makes the cross-sample refusal unreachable, so the
documents should record the widening rather than the code narrow to match — a spec edit, which this
plan's own rules keep out of an implementation step.

### Sixteen mutations, all run, all restored

Each applied from a backup with its match count asserted, tested on the affected modules, then
restored with the restore proved by `diff`.

| mutation | outcome |
|---|---|
| the reference checked after the catalog is opened | 1 test fails |
| the cohort assembled before the settings are judged (the pre-C5 order) | 1 fails |
| the reference check looks at no psp | 1 fails |
| the reference check looks at the first psp only | **survives**, and should — see below |
| the reference refusal offers no fix | 1 fails |
| the catalog compared with itself | 1 fails |
| the catalog difference always named as the scan weights | 1 fails |
| the catalog refusal offers no fix | 1 fails |
| every setting named as the seed | 4 fail |
| only the first census judged | 3 fail |
| the settings compared the other way round | **survives**, and should — the comparison is symmetric in the name it returns, so `run` and `recorded` are documentation and not behaviour |
| the pairing assertion removed | 1 fails |
| fresh samples listed rather than counted | 14 fail |
| the report ends without the command | 1 fails |
| the run's depth ladder is not the census's | 7 fail |
| the kept positions digested in reverse | 7 fail |

**Why the first survivor is not a gap.** Every psp of an opened cohort was walked against one
assembly, because the opener requires them to agree on the repeat catalog and a catalog's header
carries the whole-reference digest and the contig table it was built on. So the first psp answers
for the cohort; the loop over all of them costs one comparison a sample and stops that argument
having to hold. The method's doc says so.

### What was measured

- **The parameters file fitted over the walked fixture is byte-identical to the one from before
  C1** — 20,741 bytes, `cmp` clean against `tmp/c1_before.toml` — so nothing about a cohort that is
  not stale changed. Measured with a temporary test that wrote `file.to_toml()` to a scratch path,
  run, then deleted.
- **The affected modules' tests: 213 pass**, six of them new (three at the command, three on the
  judgement itself), one rewritten, and C4's command-level test replaced by the two that now refuse
  the same cohorts earlier.
- Gates on the final tree are in the commit message, compared as sets against the milestone's
  baseline.

**One measurement that did not happen, and how it was caught.** The first gate run of this step
produced nothing: `nohup tmp/b/gate.sh c5` failed with *permission denied*, because the script is
not executable and every previous run had invoked it through `bash`. The launcher exited 0 and its
notification said the command had completed, so the only sign was one line in the log. **A gate is
worth nothing unless the log is read**; the numbers in the commit message come from a second run,
launched with `bash`.
