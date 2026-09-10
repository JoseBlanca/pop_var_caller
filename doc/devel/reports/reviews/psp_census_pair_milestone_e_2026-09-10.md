# Code Review: the census lives inside the psp — Milestone E

**Date:** 2026-09-10
**Branch:** `census-vs-psp-perf`
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestone E
**Implementation report:** [ng_psp_census_pair_milestone_e_2026-09-10.md](../implementations/ng_psp_census_pair_milestone_e_2026-09-10.md)
**Milestone D's review:** [psp_census_pair_milestone_d_2026-09-10.md](psp_census_pair_milestone_d_2026-09-10.md)

---

## How this review is run

**Milestone A's arrangement, unchanged.** One read-only agent per step over the grouped categories,
forbidden to edit any file, to write scratch files, or to run `cargo` or any of the scripts under
review, and asked to **name** the mutations it wants as exact old→new text rather than run them.
The orchestrator runs them one at a time from a backup and proves each restore with `diff`.

**What is different at this milestone is what is being reviewed.** E1 changes no library code. Its
subjects are three shell scripts, one new example, and prose in a fourth file — so "does this
compile" is not the question, and "can this check fail" is. The brief said so, and pointed the
agent at Milestone A's own finding of a check that could not fail in the same script.

---

## E1 — the scripts, and the whole-file oracle on real reads

**Reviewed against:** the working tree over `9fdb6d6a`, four files. **One blocker.**

### Blocker — the run's headline question was asked and nothing was asserted about the answer

The script's own header says it exists to answer *what do the fitted numbers change?*, and the
`awk` block that answers it printed four sentences and set no exit status. Two regressions passed
it silently:

- **a fitted calling pass that fell back to the compiled-in defaults.** Both VCFs are then the same
  file, so the block prints `0 defaults-only, 0 fitted-only`, `0 differ`, `0 genotypes differ`, and
  the script exits 0. On this cohort the true values are 196, 3, 87 and 113, so every counter's
  correct value is non-zero and every one of them reads as zero.
- **both passes emitting a header-only VCF.** `grep -cv '^#'` prints 0 for each, inside a `$( )`
  in an `echo` whose exit status is discarded, and the block prints zeros again.

**Two checks went in, and the split between them is the point.** The `awk` block now stops the run
if either pass emitted no record at all. It deliberately does **not** demand that the two VCFs
differ: at one sample, or on ground where the fit lands close to the defaults, calling twice and
getting the same file is a legitimate answer, and a script that required a difference would fail a
correct run — the range rule in `CLAUDE.md`. What closes the first regression instead is a check on
the invocation: each pass's report line says *numbers behind the calls: N of 7 groups the file says
were fitted*, and the defaults pass must say 0 while the fitted pass must say more than 0.

A cohort too small to fit any group now stops here. That is right for this script and is said in
its comment: the script exists to answer what fitting changes, and a run with nothing fitted cannot
answer it.

### Should-fix, all applied

- **The repair's summary line was never printed.** `tail -1` was copied from the two commands
  above, whose summary is their last line. `CensusReport::lines` puts the summary **first**, and
  the log holds both streams — every sample's line goes to stderr as it is decided and again to
  stdout in the report — so in the merged log the summary sits in the middle and neither `head -1`
  nor `tail -1` finds it. It is now found by name, and a run that prints no such line fails.
- **`cp`'s exit status was unchecked**, and every downstream consequence pointed away from the
  cause: a copy that never landed makes `cmp -s` return 2, which the `if` reads as *differs*, so
  the arming check passes; the final comparison then blames the rebuild for a file that was never
  copied. `cp` is checked and the copies are counted against the psps.
- **The parameters file was never looked at.** An empty one would first be noticed by the second
  calling pass, which would report itself as the failure. It is checked for content immediately
  after the fit.
- **A test name that names nothing.** `the_two_producers_agree_on_a_cohort_with_a_repeat_tract` is
  spelled in `examples/ng_census_route_cost.rs` and in Milestone A's report and nowhere else in the
  tree. What exists is the module `ng::run::census_from_psp::the_two_producers_agree` and its three
  tests. The example now names the module.
- **"both ship" oversold the second route.** The route-cost harness's second arm walks a psp with
  an empty trailer and rebuilds from it; **no command produces a psp in that state** —
  `generate-psps` always builds a census and spec §3.3 rules out a flag to skip it. The rebuild
  half is real, and it is the repair's pass. Both the harness and its wrapper now say that: one
  route is what every run does, the other is what repairing a psp that arrived without a census
  costs, and the harness makes that state deliberately in order to price it.
- **`ng_census_agreement_mutations.sh` reported its own failures as survivals.** Under `set -u`
  alone, its `cp` and its python replacement were both unchecked, and the python contains
  `assert find in s`. A drifted anchor therefore left the file unmutated, the tests then passed,
  and the harness printed `test result: ok` — a mutation harness saying *nothing catches this*
  about a defect it never wrote. Both statuses are checked and a failure stops the run.
- **And its fourth case's expectation was stale.** It said a change to a read's minted error cannot
  move a census byte, because a census holds depth codes and allele counts and no per-read quality.
  That stopped being true when the fit stage's Milestone C put the per-read-group minted read-error
  totals into the census. Measured rather than reasoned: applied on its own, the fourth defect
  fails both of the module's comparing tests. The reason the third test does not fail is not that
  its census has no read in it — it is built on the same sample as the first — but that it never
  compares the two censuses; it checks the fixture's walk reached its kept loci.

### Minor, applied

- **Both `cmp` loops wrote their own pair of paths.** One function now answers *does the walked psp
  differ from its copy*, and both loops call it, so the second cannot drift into comparing a copy
  with itself. The mutation that does exactly that is in the suite below.
- **The walk's sample count is checked against the psps on disk**, so the heading *the walk left
  one file a sample* is now a claim the script tests rather than one it asserts.
- **Each number parsed out of the walk's report must match exactly once.** The per-sample lines say
  "of which N bytes are **its** census" and are printed twice a sample, so a pattern one word
  looser returns twelve numbers for six samples — the shape this check was already found in once,
  at Checkpoint A. Twelve numbers would make `[[ -z ]]` false and `(( inside == 0 ))` a bash syntax
  error evaluating false, and the run would carry on.
- **`ng_psp_drop_census`'s failure message re-invented a worse one.** `TrailerReplacementFailure`
  already says whether the file is as it was or cut short and footerless — the difference between
  running the tool again and re-walking the sample — and the message printed the enum's `Debug`
  instead. It now prints the failure's own words and its cause, and names the psps it had already
  emptied.
- **One sentence of that example's doc was wrong about the format**: it said the footer's trailer
  offset *and* length move with the truncation. Only the length moves.
- **An ambiguous antecedent** in the route-cost harness: *the second route* meant the harness's own
  arm everywhere else in the file and the shipped command in that sentence. The command is named.

### Nits not taken

- **`rm -rf "$out"` on an unvalidated fourth argument**, at the top of
  `ng_fit_stage_end_to_end.sh`. Pre-existing, not part of this change, and left alone.

### What the reviewer confirmed rather than faulted

- **The oracle is armed, and there is no path to the final comparison without a rebuild.** All
  three ways it could be hollow — copies handed over unemptied, the originals emptied by mistake,
  the repair skipping instead of rebuilding — leave one of the two comparisons failing.
- **No shipped route writes a census with a pileup identity.** There are exactly two `write_census`
  call sites outside tests and examples, the walk's and the repair's, and both pass `None`.
- **The route-cost harness's in-clock census write is matched by real work.** `regenerate-census`
  encodes the census and hands it to `replace_trailer`, which writes the same bytes back into the
  psp — one write of a census either way.
- **The new example belongs in `examples/`, and no smaller route exists.** `dd` cannot do it,
  because the footer has to be re-encoded behind the cut; `--force` on the repair is ruled out by
  spec §8; a `#[cfg(test)]` helper is unreachable from a shell script.

### The one limit that stands, and it is inherent

**The arming check is the only thing between this oracle and a vacuous pass.** If the copies are
never emptied *and* the check that notices is also disabled, the repair skips every psp, the copies
stay byte-identical to the originals, and the final comparison passes. The mutation suite below
runs that pair deliberately and records it as surviving. It is not a defect to fix in the script:
the before-comparison *is* the guard, and a guard cannot guard itself.

---

## E2 — the sidecar's machinery deleted

**Reviewed against:** the working tree over `fe7431f9`, 16 files, about 413 lines deleted and 154
added. **No blocker.** The reviewer's own summary: the deletion is correct and complete, and every
finding is in prose that survived it, plus one untested branch.

### What the review confirmed rather than faulted

- **Nothing deleted was doing work.** Checked against the pre-change tree: every reader of
  `CensusFile.pileup` was a test assertion, and `freshness`, `freshness_by_header` and `Freshness`
  had no reader outside their own module at all.
- **The format is safe.** The flag byte is the header's last field, so a census written with the
  identity absent — every trailer since Milestone A — encodes to the same bytes and decodes to the
  same value. `VERSION` is right to stay at 4.
- **`psp_beside`, `CensusInCohort` and `open_census_cohort` were already gone**, deleted at C2 and
  named in that commit's own message. E2's list was partly satisfied before the step began.
- **A small saving nobody had claimed:** dropping the header digest from `census_from_psp` removes
  a second file open and full header read per psp, on top of the reader open that follows.

### Should-fix — the untested branch, and the fix changed the behaviour

The arm that stepped over a present identity **could not be reached by any writer in the tree**,
and the test that had covered the format's two flag values went with the type. Both of the
reviewer's routes were open: test it with a hand-built census, or refuse. **Refusing was taken**,
for the reviewer's own argument — such a census names a psp by values this build no longer compares
against anything, so reading past them accepts a file it can say nothing true about — and a test
now builds one, finding the flag byte from `encode_header`'s own output rather than guessing an
offset.

### Should-fix — the version had nothing behind it

Every test compares a file's version word against `VERSION`, so all of them pass whatever it holds,
and this step's central claim — that the version did not move — had no test. A second constant and
a test that the two agree makes a bump a deliberate edit in two places. The reviewer's broader
point, that there is **no golden census on disk**, stands and is bigger than this step: a
coordinated change to both sides of the flag byte would pass every test while making every psp
already written unreadable. It goes to Checkpoint E.

### Should-fix — nine sentences that had become false, all applied

Six in the files this step touched, each promising in the present tense a check the step deleted:
`gatherer.rs`'s two, `generate_psps/tests.rs`, `regenerate_census/tests.rs`,
`regenerate_census.rs`, and `census_from_psp.rs`'s digest test. Three outside them, each a live
engineering rationale resting on the deleted mechanism: the psp header's read-filter keys, the
writer line's "encoded twice for now", and the route-cost harness's constant command line. All
rewritten around the reason that still holds — for the filters, that a cohort's psps are refused
unless their headers agree.

**One test was renamed rather than repaired.** `the_census_it_writes_into_the_trailer_names_no_psp`
had been reduced to a bare decode; it is now
`the_trailer_it_writes_decodes_as_a_census`, which is what it checks.

### Minor, applied

The struct's summary line still said "which psp it came from"; the module's format table had lost
the read groups, the minted totals and the flag byte; and a comment explaining a discard sat above
a line that no longer discarded anything.

### Two documents outside this plan, left alone

`doc/devel/ng/arch/run_streaming.md` cites `PileupIdentity::of_header` and `freshness` by line, and
`doc/devel/ng/spec/run_streaming.md` §6.1 justifies what the psp header must contain by a consumer
that no longer exists. Both are present-tense claims about the tree and both are now false. Editing
another plan's spec is not this plan's to do — they are listed at Checkpoint E.

### What is not guarded any more, and was not before either

The rule that a census's staleness must never key on a modification time (spec §6.1) lost its only
test with the two functions. The reviewer checked that nothing in `census_freshness.rs`,
`regenerate_census.rs` or `psp_caller.rs` reads a timestamp, so the property holds structurally —
but it holds by nobody having written the code, not by a test. Noted at Checkpoint E.
