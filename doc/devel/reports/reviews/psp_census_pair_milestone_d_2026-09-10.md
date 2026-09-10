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
