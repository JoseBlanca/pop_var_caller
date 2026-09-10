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

---

## D2 — fresh psps skipped

**Committed:** see `git log` for `feat(ng): D2`.

### What it does

**A psp whose census is the one this run would write is left alone, and the run says so.** What
counts as needing nothing is every reason a census cannot be used — spec §4.2's three causes, and
damage of the trailer with them — because this command reads the reference and rebuilds the
selection anyway, so comparing each census against them costs one open and one read a psp rather
than a second pass over the genome.

So a run stopped part-way and started again does only the samples still owed, and a build whose
selection constants changed rebuilds every one without being told to.

**The order the run now takes**, which is also spec §8's and the finding carried out of D1's review:

1. the `--psp` arguments expanded;
2. **the cohort opened** — its agreement on the ground, the catalog and the criteria is what
   everything else rests on, and a mistyped path is answered here rather than after a reference
   read;
3. **the psps' heads judged**, before the reference is read;
4. the reference read, and compared with every psp header;
5. the catalog opened and compared with the one the psps record;
6. the segments cut and the selection planned;
7. **each fresh head's census compared with what this run records under** — the third cause;
8. the cohort closed, and only the psps owed a rebuild rewritten.

**What the head pass buys, since the census read reaches the same two causes on its own.** An empty
trailer and a census of another format both come back from the census read — it checks the same
magic and the same version word — but only after the reader's head read has taken up to a mebibyte
to get there, where the head answers in ten bytes. It changes no outcome; it changes what a cohort
of stale psps costs to judge. The step's review measured that: with the head's verdict forced to
*fresh*, the whole suite still passes.

**Two judgements taken deliberately, both recorded in the code.** A psp that will not read is owed
a rebuild rather than refused here — elsewhere the two are kept apart, because regenerating the
census of an unreadable file fixes nothing, and here the distinction dissolves, since this
command's own work is to read that psp and the attempt names it if it cannot. And a census whose
head is this build's but whose sections will not decode is owed a rebuild too: no cheap read can
tell it from a whole one, and rebuilding rewrites exactly the bytes that are damaged.

### What it costs, measured rather than estimated

**A psp that needs nothing costs one open and up to a mebibyte read, of which a few hundred bytes
are decoded** — the census reader's head read, and less than that for a census shorter than it.
Against a whole psp's records, which is the pass this command exists to avoid, it is a rounding
error; over a thousand samples it is about a gibibyte, which is worth knowing. **The first draft of
this section said "a few hundred bytes" for the read itself**, which is what gets *decoded* — wrong
by about 3,000× — and the review caught it in three places.

### How "its records are not read" is shown

One psp's blocks are corrupted before the run: 32 bytes overwritten just below where the index
begins, which is inside the last block. A run that read that psp's records would fail; this one
succeeds and skips it.

**The control is in the same test**, and it is what makes the first half mean anything: with that
psp's census emptied, the same corrupted file is owed a rebuild, and then the run does fail on it,
naming the sample and the file. Without the control, a psp that would have read cleanly anyway
would pass the first half.

### What D2 exposed in D1's tests

**Over a freshly walked cohort, every test in the file now compares the walk's own bytes with
themselves** — because a freshly walked cohort is exactly what this step skips. Seven of them empty
the psps' trailers first, which is D1's blocker fix generalised: the bytes a comparison asserts have
to be bytes the run under test wrote.

### What the review found, and what was done

**No blockers.** Beyond the wrong size above: the error type's doc still described the order this
step replaced and still said D2 would make the re-run cheap, in the commit that makes it so; the
report's "put nothing into the fit" line divided by the rebuilt count without saying so; the skip
list needed the argument for why it lists what the refusal report only counts; and four tests could
not fail — the report's singulars and its no-skips line were never rendered, the selection test used
the budget where the kept positions differ anyway, and nothing pinned that the cohort is opened
before the reference is read. All fixed, the last with a test where a bad `--psp` and a missing
`--reference` are wrong at once, so whichever is read first is the one that refuses.

**One finding routed to a later step.** Plan step D4's whole-file oracle — a walked psp copied,
regenerated, and compared byte for byte — now passes without a rebuild happening, because the copy
is skipped. D4's own line in the plan records that its copy's trailer has to be emptied first.

### Thirteen mutations, eleven measurements

| mutation | outcome |
|---|---|
| the head's verdict is always fresh | **survives**, and that is the finding: the census read reaches the same causes, so the head pass is about cost |
| a census that will not decode is treated as fresh | 1 test fails |
| a psp whose head will not read is treated as fresh | **survives** — nothing in the tree can produce a head-read failure |
| the recorded settings are not compared | 1 fails |
| only the kept positions are compared | 1 fails — the half-budget arm the review asked for |
| every stale head is treated as fresh | 9 fail |
| the verdicts are paired with the wrong psps | 11 fail |
| a run that skipped nothing says it skipped none | 1 fails |
| the skipped count is always plural | 1 fails |
| the rebuilt count is always plural | 1 fails |
| the corruption lands on the block index | 1 fails — which is what says the 32 bytes are in the blocks and not in anything the open pass reads |
| the reference is checked after the catalog is opened | 1 fails |
| the cohort is opened after the reference is read | 1 fails — the order carried out of D1 |

**Two of the thirteen were faulty on the first run and are not counted as measurements**: one did
not compile, and one moved the cohort block back to where it already was, so it proved nothing. Both
were fixed and re-run, and both are caught. **A mutation that does not compile, or that changes
nothing, looks exactly like one no test catches** — the driver prints the apply step and the test
result for each, which is what made the difference visible.

### And one report of mine that was wrong

Before those two were re-run I told the owner the thirteen were running when the driver had died on
a syntax error I had introduced in it — the mutation list's closing quote lost its newline, bash
rejected the file, and my check was a peek at an empty log, which cannot tell *not started* from
*just started*. The relaunch was verified three ways: the log grew, the first mutation reported its
substitution applied, and the driver process was alive.

---

## D3 and D4 — a stopped run, and the whole-file oracle

**Committed:** see `git log` for `feat(ng): D3+D4`. **One loop iteration and one commit for both
steps**, which the plan-driven skill allows for tightly-coupled adjacent steps and which is named
here: both are test-only, over the same command, and both rest on the same idiom — a psp whose
census this run would write is skipped, so a test that wants a rebuild has to make the psp owed
first.

### D3 — a run stopped part-way does only what is left

D2's skip rule is what makes this true; the step is the test that says so. The first run rebuilds
`alpha` and then fails on `zeta`, leaving the cohort half repaired. The bytes are put back, and the
second run skips `alpha` and rebuilds `zeta` alone.

**The property nothing else covers: the second run skips a census *this command* wrote.** Every
other skip test skips the bytes the walk wrote. If what `write_census` records and what
`CensusPlan::recording_terms` is compared against came apart, this command would rebuild its own
output for ever and a stopped run would never converge, however many times a person ran it.

**How the psp is made to fail, and why not the way the plan asked.** Its blocks are corrupted — 32
bytes just below where the index begins — because the rebuild is the only pass that reads records.
Permissions would be plainer and cannot be used: this suite runs as root inside the dev container,
where a read-only file is not read-only. The corruption stands in for a read that failed and then
stopped failing, which is why the bytes go back between the runs.

**Two psps rather than the plan's three**, recorded in the plan's own D3 line: the fixture cohorts
have two samples, and the property needs one of each. **What two cannot see** is a run that pressed
on past the failure and rebuilt *later* samples before returning the first error — with the failing
psp last there is nothing after it to have been touched. Closing that needs a three-sample fixture.

### D4 — a regenerated psp is the walked one, byte for byte, whole

Both psps of the fixture with a repeat tract are copied into a directory of their own, their
trailers emptied, and the command run over them; each file then equals the walked one byte for byte
— header, blocks, index, trailer and footer.

**The emptied trailer is the plan's own note on this step**, added at D2: a psp whose census this
run would write is skipped, so a copy handed straight to the command comes back untouched and a byte
comparison passes without a rebuild having happened. It also makes the comparison stronger than the
plan asked — the footer's trailer length goes to zero and has to come back.

**What it covers that this file's trailer comparisons do not**, and it is two things rather than the
four the first draft claimed: the header's exact **bytes**, where the sibling test compares its
`Debug` rendering, and every record's **content**, where that test compares only how many there are.
The review corrected the other two: a dropped block is caught by that record count, and a rewritten
index never reaches any comparison because the reader refuses the file first.

**And the honest limit, which is now in the test's own doc.** Against this command as it stands there
is no small defect this catches first — the only write is `replace_trailer`, which touches nothing
below the trailer's offset. What it is for is the change that would stop that being true: a rebuild
that wrote the psp through `PspWriter` again, re-compressed its blocks, or re-encoded its header.

### What the review found, and what was done

**No blockers**, and it confirmed neither test can pass on an idle command. Fixed: the over-claimed
coverage above; "three read groups" said of a psp that declares two (the *cohort* has three, which
is why both psps are now copied); the missing guard that the comparison covers the census's
stratum-keyed half at all, now a tally assertion; the reason permissions were not used; the claim
that nothing re-reads a skipped psp's records, which this test cannot see and D2's can; a quarter of
an hour without its subject; and a variable called `whole` that held a psp with an emptied trailer.
The corruption, written out twice, is one helper now.

### What was measured

- **6,738 lib tests pass**, against 6,736 at D2: these two steps are two tests.
- Mutations below. Gates in the commit message, as sets against the milestone baseline.

### Ten mutations, nine caught — and the survivor corrected a claim of D2's

| mutation | outcome |
|---|---|
| the rebuild order is reversed | 2 tests fail, and only the stopped-run test could see it |
| a sample that will not rebuild is passed over instead of stopping the run | 2 fail |
| the tail is never rewritten | 9 fail |
| nothing is ever skipped | 4 fail |
| everything is skipped | 13 fail |
| the census names the psp it came from | 8 fail |
| only the first psp owed is rebuilt | 9 fail — the reason D4 copies both psps rather than one |
| **a skipped psp's records are read anyway** | **survives** |
| the fixture's repeat tract falls below the period-2 floor | 3 fail, one of them the tally assertion this step added |
| the tail rewrite drops the index checksum | 13 fail |

**The survivor is a correction, not a gap to close.** D2's report and this command's module doc said
a psp needing nothing "keeps its records unread". What the corrupted-block test pins is narrower and
is the part that matters: **no rebuild pass happens for a skipped psp** — a rebuild reads every
record and would fail on the corrupted ones — and that pass is the quarter of an hour a sample.
Code added to the skip arm that reads the records and discards the result passes every test in the
file, because a discarded failure is invisible and nothing counts block bytes the way
`trailer_bytes_read` counts trailer bytes. Both the module doc and the test's own doc now say the
narrower thing, and say how it was measured.

**Two of the ten aim at code this milestone did not write** — the fixture's tract length and the
trailer rewrite's footer — and they are here because D4's oracle is the only test in the tree that
would notice either. Both are caught, the second by thirteen tests, which is the evidence for the
review's point that the reader refuses a disturbed footer before any comparison reaches it.
