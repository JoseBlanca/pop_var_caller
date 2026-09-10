# Code Review: the census lives inside the psp — Milestone D

**Date:** 2026-09-10
**Branch:** `census-vs-psp-perf`
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone D
**Implementation report:** [ng_psp_census_pair_milestone_d_2026-09-10.md](../implementations/ng_psp_census_pair_milestone_d_2026-09-10.md)

---

## How this review is run

**Milestone A's arrangement, unchanged through C.** One read-only agent over the grouped categories,
forbidden to edit any file, to write scratch files, or to run `cargo`, and asked to **name** the
mutations it wants rather than run them. The orchestrator runs them one at a time, restoring from a
backup and proving the restore with `diff`.

**One thing went differently and it is worth recording.** D1 was committed while the agent was
still reading, on the owner's instruction, so the review arrived against a commit rather than a
working tree — and its findings land in their own commit on top. The agent noticed and re-cited
everything to the committed text. It also reviewed the two files the brief had not listed, the
implementation report and the plan's tick, and two of its findings are about the report.

---

## D1 — the rename, and the write into the psp

**Reviewed against:** `34359259`, twelve files. **One blocker.**

### Blocker

**B1 — nothing in the suite could fail if the command never wrote anything.** `replace_trailer` is
the whole of what D1 adds, and every test that could have observed it compared a psp's trailer with
the trailer *the walk had already put there*: the two parity tests captured the trailers before the
run and asserted they were unchanged after, the identity test decoded a trailer that was the walk's
either way, the cohort-assembly test read the walk's censuses back, and the two refusal tests
asserted that no psp was rewritten — which also holds when no psp is ever rewritten at all. **The
review measured it: with the write removed, all fourteen tests passed.** Plan step D4's whole-file
comparison would have passed on that same no-op, so the milestone as planned never proved the write
lands.

*Fixed*: one test, `the_census_is_written_into_a_psp_that_has_none`. It empties one psp's trailer
before the run and asserts that psp comes back carrying the census its walk had written, while the
other psp is left as the walk sealed it. The technique was already in the tree —
`estimate-parameters`' own tests empty a trailer with `replace_trailer(path, b"")` to make a psp
carry no census. The mutation table below is the measurement that it bites.

### Should-fix, all fixed

**S1 — spec §8 asks for the cohort's agreement to be checked first, and it is checked fourth.** The
run reads the reference and turns its unambiguous runs into selectable regions before it opens the
cohort, so a mistyped `--psp` path costs a full reference read first, and the same person meets the
two commands' refusals in two different orders — `estimate-parameters` opens its cohort first.
*Not fixed here, and it is the one finding carried:* the reorder is a change to the shape of the
run rather than to a message, the review's own mutation 12 predicts it costs no test, and D2 has to
rework this preamble anyway to judge freshness before rebuilding. **Carried to D2**, where the skip
rule lands in the same block.

**S2 — the error type's doc claimed a pre-flight the command does not have.** It said every refusal
but one comes before a psp is rewritten. Three arrive inside the loop, and for sample *k* they
arrive after *k*−1 psps have been rewritten. The command this replaced genuinely had that property
and could: it judged every output path before doing any work, because what it wrote was a separate
file. *Fixed*: the doc now says a cohort of sixty that fails at the fortieth leaves thirty-nine
rebuilt, that no report is printed on failure so the per-sample progress lines are the record of
what was done, and that plan step D2 is what makes the re-run cheap.

**S3 — a bare intra-doc link to a type the module does not import.** `broken_intra_doc_links` is
denied in `Cargo.toml`, so `cargo doc` would have rejected it. *Fixed* with the full path.

**S4 — the help text offered a cause no shipped command can produce.** It listed "one whose census
an append discarded"; spec §3.4 records the owner's ruling that no shipped command calls `append`
and there is no user case for it. The text it replaced named the case that is real — a psp written
before the census moved into the trailer, spec §4.2's *no census* row. *Fixed* in both the help and
the module doc.

**S5 — three doc comments in `src/` that this commit falsified**, each true at `93fb1ac3`: that
`generate-census` still writes census files and `estimate-parameters` still reads them
(`generate_psps.rs`); that `census_from_psp` computes an identity "because `generate-census` still
writes census files that need one" (`census_from_psp.rs`); and the list of commands that refuse a
disagreeing cohort (`run/mod.rs`). *All three fixed*, and the second now says what is true: the only
part of that identity any caller reads is the record count, and the header digest beside it is a
second open-and-hash with nothing left to compare against.

**S6 — the codec test's stated caller no longer does what the sentence said.** Nothing in the tree
decodes a census and re-encodes it any more. *Fixed*: the test keeps its place and its reason is now
that it guards the codec against a lossy decode landing unnoticed, which no round-trip test over
decoded *values* can see.

**S7 — the record count the report prints was pinned against nothing.** It came from the psp and the
only test that read it built the expected string out of the same field, so `records: 0` passed —
which the deleted test had said in as many words. *Fixed*: the line's count is asserted against the
psp's own record count, read from the file.

**S8 — two tautologies and an assertion that passes when clap does *not* answer to the name.**
`SUBCOMMAND` is defined as the constant, so `assert_eq!(SUBCOMMAND, THE_COMMAND_THAT_REBUILDS_A_CENSUS)`
is `X == X`; and the `--help` check asserted only that the error text contains the word, which
clap's *unrecognised subcommand* error also echoes back. *Fixed*: the tautologies are gone with a
comment saying why they were, and the `--help` check asserts `ErrorKind::DisplayHelp`.

**S9, S10 — two claims in the implementation report.** The parity section's "what it no longer
covers" named two things that are covered elsewhere and missed the two that were genuinely
uncovered (the write itself, and the record count); and "nothing in the tree still names the command
that was renamed" was false in the same paragraph that went on to name three things that do.
*Both fixed in the report.*

### Minor, fixed

`run_ground`'s two doc comments named one psp-taking command where there are now two, and the
defensive-arm argument holds for both. The order sentence in the error doc listed the selection
last where one of its own refusals fires second. The module doc now says what an interrupted repair
costs — the tail is truncated before the new census is written, so a psp caught between those two
moments has no footer, cannot be read, and therefore cannot have its census rebuilt from its
records: that one sample has to be re-walked from its alignments. The command this replaced could
only destroy a cache file beside the psp. The reference's `Arc` is gone (never cloned, never crosses
a thread) and `Segmentation` is imported like every other `ng::run` type in the file.

### Minor, recorded and not fixed

- **An incomplete psp — a walk killed part-way, no footer — is reported as `Unchanged`**, whose
  message says the call can simply be made again, and it will fail identically for ever. The cause
  underneath says the file is incomplete, so the information is in the chain; only the instruction
  misleads. It belongs to `replace_trailer`'s own verdict rather than to this command.
- **The seven-step preamble is duplicated** between this command and `estimate-parameters`, about
  55 lines including two argument-identical catalog-check wrappers. Both are pinned to the walk by
  separate oracles, so a divergence would be caught rather than silent. Worth a note for E2.
- **`CENSUS_FILE_EXTENSION` and `census_path_for` now have no caller in `src/`**, which is exactly
  the state plan step E2 waits for.
- **`scripts/ng_fit_stage_end_to_end.sh` is now broken**, not merely out of date: it invokes
  `generate-census` with `--output-dir`. Step E1 owns it and depends on D4.

### Nits, fixed

`expect_err("{flag} …")` is not a format macro, so the braces printed literally — and fixing it
found a type error in my own first attempt at the fix. "Its second sample has no reads" was true in
the fixture's alignment order and false in the order this command reads them, which is name order,
so it now names the sample. "A quarter of an hour a sample" appears with its subject: 0.25 to 3.25 s
on the fixtures spec §2 measures, a quarter of an hour at 50× human. The plan's D1 line said
"fifteen tests" where the module had thirteen.

### What the review said about the tests, and it is the useful half

**Strongest:** the sample-to-psp pairing, the plan's identity (the parity tests do bite on a
selection divergence — mutation 8 below), both file refusals, and the two fixture-shape assertions
that stop the parity tests passing over two empty halves. **Weakest, and now fixed:** everything
that rested on a trailer comparison alone (B1), the record count (S7), and the name assertions (S8).

### The mutations

Nine, each applied from a backup with its match count asserted, tested on the affected modules, and
restored with the restore proved by `diff`. **All nine were caught**, which they were not before the
fixes above: the first is B1's own measurement.

| mutation | outcome |
|---|---|
| the write never happens | 1 test fails — the new one, where all fourteen passed before it |
| the census keeps its pileup identity | 4 fail |
| each psp is paired with the wrong sample | 2 fail |
| the reference is not checked against the psps | 1 fails |
| the catalog is not checked against the psps | 1 fails |
| the criteria outrank the reference in the catalog comparison | 1 fails |
| the reported record count is a constant | 1 fails — the assertion S7 added |
| the selection uses another seed | 3 fail, including both parity tests |
| the library advertises the old command name | 1 fails, and the binary aborts: one test parses the printed word with `Cli::parse_from`, which exits the process on a name clap does not know |

**Two of the review's own mutations were not run, and both were arguments rather than measurements.**
Comparing the two catalog headers the other way round is symmetric in five of its six clauses, so it
would survive and prove nothing; and holding the cohort open across the writes cannot be observed
from a test on Unix, where an open descriptor does not stop a truncation. Both are recorded here
instead: the second means *closed before the first write* is a design property with no test behind
it, and spec §5's descriptor and memory ground is what argues for it.

---

## D2 — fresh psps skipped

**Reviewed against:** the working tree over `8f3485eb`, four files. One read-only agent over the
grouped categories. **No blockers**, and it confirmed the three things the brief asked it to
challenge: the skip rule reaches all of §4.2's causes and maps each to the right action, the new
order keeps the reference check ahead of the catalog open, and *its records are not read* is true —
nothing between opening the cohort and the decision touches a block.

### Should-fix, all fixed

**S1 — the head judgement changes no outcome.** An empty trailer and a census of another format
both come back from the census read as well, since it checks the same magic and the same version
word, and both arrive as an error that already means *owed a rebuild*. The review's own mutation
confirmed it: with the head's verdict forced to *fresh*, the whole suite passes. *Kept, with the
reason rewritten*: what the head buys is that those two causes are reached in ten bytes rather than
after a read of up to a mebibyte. It changes what a cohort of stale psps costs to judge, not what
the run decides — and the comment says exactly that now.

**S2 — a size claim wrong by about 3,000×, in three places.** "A few hundred bytes" is what gets
*decoded* out of a census's head; the read itself is `HEAD_READ_BYTES`, a mebibyte, so about a
gibibyte over a thousand samples. The file's own older sentence had it right — "at most a megabyte,
and less for a census shorter than that — out of which the census's header and its directory are
decoded, a few hundred bytes of them" — and the new comments had kept the second half and dropped
the first. *Fixed*, and each now says both numbers and that it is a rounding error against the
record pass it avoids.

**S3 — the error type's doc described the order this step replaced**, and still said plan step D2
was what would make the re-run cheap, in the commit that makes it so. *Fixed.*

**S4 — both `Err(_) => true` arms swallow an i/o failure**, where the module that produced the
verdict draws the line the other way: a psp that cannot be read is not a stale psp, because
regenerating it fixes nothing. *Recorded rather than branched on*: here that distinction dissolves,
because this command's own work is to read the psp, so an i/o fault reaches the rebuild, which reads
the file and names it. The comment says so, for both arms.

**S5 — a test claimed to be the only one that observes the write**, which stopped being true when
seven others began emptying their trailers. *Reduced to what it still covers*: the pair — one psp
rebuilt and one skipped in the same run.

**S6 — the selection test used the easy budget.** At three positions the kept set differs too, so it
passed against a check comparing the kept positions alone; the review's mutation proved it. *Fixed*
by looping over half the shipped budget and three, which is what the sibling test at the fit does
and for the same reason.

**S7 — the report's singulars and its no-skips line were never rendered.** *Fixed*: the mixed run
asserts "regenerated 1 census … and skipped 1 psp that needed nothing", and the all-rebuilt run
asserts the first line says nothing about skipping.

**S8 — the milestone bookkeeping.** *Done in the commit*: the plan's tick and both report sections.

**S9 — plan step D4's oracle is now a no-op**, because the copy it compares is skipped. *Recorded in
the plan's own D4 line*: the copy's trailer has to be emptied first.

**S10 — a doc sentence with no verb.** *Fixed.*

### Minor, fixed

The "put nothing into the fit" line divides by the rebuilt count, and now says so. The skip list
needed the argument for why this report lists what the refusal report only counts — that one is
about what a person must go and do, this one is the record of what a run did — including that a
re-run of sixty samples prints sixty lines. "All three of §4.2's causes" left out trailer damage,
which the code also handles. "The third cause costs nothing" was wrong about what is free: the
reference read and the selection rebuild are, the comparison is one open and one read a psp. The
verdicts and the readers are now paired under a `debug_assert`, and a trailer holding another
sample's census is owed a rebuild — one string comparison, and the only thing that would catch a
spliced file.

### Recorded, not fixed

- **Closing the cohort before the writes has no test behind it**, and cannot have one on Unix, where
  an open descriptor does not stop a truncation. Its argument is memory: a thousand descriptors and
  a thousand block indexes held while the psps are rewritten (spec §5).
- **The progress stream and the report order their lines differently** — skips are printed as they
  are decided, before the first rebuild, and the report lists rebuilds first.

### The mutations

Thirteen, each applied from a backup with its match count asserted and restored with the restore
proved by `diff`. **Eleven caught, two survivors, both predicted.**

| mutation | outcome |
|---|---|
| the head's verdict is always fresh | **survives** — S1's own measurement |
| a census that will not decode is treated as fresh | 1 test fails |
| a psp whose head will not read is treated as fresh | **survives** — no head-read failure is producible |
| the recorded settings are not compared | 1 fails |
| only the kept positions are compared | 1 fails |
| every stale head is treated as fresh | 9 fail |
| the verdicts are paired with the wrong psps | 11 fail |
| a run that skipped nothing says it skipped none | 1 fails |
| the skipped count is always plural | 1 fails |
| the rebuilt count is always plural | 1 fails |
| the corruption lands on the block index | 1 fails |
| the reference is checked after the catalog is opened | 1 fails |
| the cohort is opened after the reference is read | 1 fails |

**Two of them were faulty on their first run and produced no measurement**: one did not compile, and
one moved a block back to where it already was. Both were fixed and re-run. **A mutation that does
not apply, does not compile, or changes nothing is indistinguishable from one nothing catches**, and
the only reason it was visible here is that the driver prints the apply step and the compiler's
output beside each test result.

---

## D3 and D4 — a stopped run, and the whole-file oracle

**Reviewed against:** the working tree over `c55c55d5`, one file — both steps are test-only. One
read-only agent. **No blockers**, and it agreed the two steps belong in one commit.

**What it confirmed, which is the half worth having.** Neither test can pass on a command that does
nothing: the stopped-run test's `expect_err` needs a rebuild to have been attempted, and its
`alpha` assertion reads a trailer the test itself emptied, so only this run's write can put the
walk's census back; the byte comparison asserts the copies differ *before* the rebuild and that the
samples were rebuilt. It also traced why the failure lands after the first psp — the arguments are
sorted, the owed list is built in that order, and the rebuild loop follows it.

### Should-fix, all fixed

**S1 — the byte comparison's unique coverage was over-listed.** It claimed a rebuild that
re-encoded the header, dropped a block, rewrote the index or disturbed the footer would pass the
trailer comparisons and fail here. Two of those are false: the record count catches a dropped block,
and a rewritten index never reaches the comparison because `PspReader::open` re-checks the index's
checksum and refuses the file. *Fixed* to the two that are true — the header's exact bytes and every
record's content — **with the honest limit beside them**: against this command as it stands nothing
small is caught here first, because the only write is `replace_trailer`; what it is for is the change
that would stop that being true.

**S2 — "three read groups" was said of a psp that declares two.** The cohort has three, split two
and one. *Fixed by copying both psps*, which gets all three, a sample with one group, and a cohort of
more than one — where a run that rebuilt only its first psp would have passed over a cohort of one.
The same sentence claimed the comparison covers the census's stratum-keyed half, which nothing
guarded; a tally assertion now does, as the sibling parity test's does.

**S3 — the plan still asked for three psps.** *Recorded in its D3 line* as an *as built* note, with
why permissions could not be used and what two psps cannot see.

**S4 — why the psp was corrupted rather than made read-only** was missing, which is what makes the
shape honest. *Fixed*, in the test and in the helper.

### Minor, fixed

The claim that nothing re-reads a skipped psp's records, which this test cannot see — it cites D2's,
which corrupts the psp it expects skipped. The strongest property said plainly: the second run skips
a census *this command* wrote. The quarter of an hour with its subject. `whole` renamed, since it
held a psp whose trailer had just been emptied. The error match pins the file as well as the sample.
`breaks` and `alone` renamed to say what they hold. A stale comment naming the oracle by its plan
step, where the test now exists. And the corruption, written out statement for statement in two
tests, is one helper.

### Recorded, not fixed

- **A run that pressed on past a failure and rebuilt later samples is invisible to a two-psp
  fixture**, since the failing psp is last. The error type's own doc claims the sixty-sample shape;
  closing it needs a three-sample fixture.
- **Nothing tests that a psp *after* the failure is left alone**, for the same reason.

### The mutations

Ten, each applied from a backup with its match count asserted and restored with the restore proved
by `diff`. **Nine caught, one survivor — and the survivor was worth more than the nine.**

| mutation | outcome |
|---|---|
| the rebuild order is reversed | 2 tests fail; only the stopped-run test could see it |
| a sample that will not rebuild is passed over | 2 fail |
| the tail is never rewritten | 9 fail |
| nothing is ever skipped | 4 fail |
| everything is skipped | 13 fail |
| the census names the psp it came from | 8 fail |
| only the first psp owed is rebuilt | 9 fail |
| **a skipped psp's records are read anyway** | **survives** |
| the fixture's repeat tract falls below the period-2 floor | 3 fail |
| the tail rewrite drops the index checksum | 13 fail |

**The survivor.** The review predicted it would survive the stopped-run test and be caught by D2's
corrupted-block test. It survives that one too: the code it adds reads the records and discards the
result, and a discarded failure is invisible — nothing counts block bytes. So the claim *a skipped
psp's records are not read* was never pinned; what is pinned is that **no rebuild pass happens**,
which is the expensive one. The module doc and the test's doc now say that, and cite the measurement.

**Two mutations aim outside this milestone's code** — the fixture's tract length and
`replace_trailer`'s footer — because D4's byte comparison is the only test that would notice either.
The second failing thirteen tests is the evidence for S1's limit: a disturbed footer is refused by
the reader long before any byte comparison sees it.
