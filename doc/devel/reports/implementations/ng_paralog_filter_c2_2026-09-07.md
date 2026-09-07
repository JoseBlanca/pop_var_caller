# The hidden-duplication filter — C2: where a finished record goes, and the flag that decides

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone C, step C2
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.4, §3.6, §1.1 goal 7
**Branch:** `ng-paralog-filter`

## The answer

**A finished record now has two places it can go, and the run with the filter off still goes
through the line it went through before this work existed — unchanged, not re-derived.**

That last word is the whole of the step. The standing oracle for the filter is that
`--paralog-fdr 0` reproduces the previous run byte for byte, and **what carries it is that the
off path still ends in the same call**: `CalledRecordSink::StraightToTheVcf` holds the
`VcfWriter` and its `accept` does nothing but `writer.write_record(record)` — no encode, no copy,
no second branch. Six lines change in each subcommand, in three edits: the writer is built *inside* the sink
rather than before it, the closure calls `accept` instead of `write_record`, and the sink closes
itself rather than handing the writer back.

**An earlier version of this step left the sink unwired**, so the two subcommands still called
`write_record` directly and the diff on the writing path was empty. That made a tidier claim and a
worse step: the sink was dead code, and C4 would have had to wire it and prove the off path all
over again. It is wired now, and the byte-identity rests on `accept`'s off arm being a single
delegation rather than on the diff being absent.

With the filter on, a record's line is parked on the spill with what the scorer will need beside
it, and the VCF is not opened at all.

## What was added

| file | what it is |
|---|---|
| [`pass_one.rs`](../../../src/ng/run/paralog_filter/pass_one.rs) | 171 lines: `entry_for`, `CalledRecordSink`, `SpillingSink`, `PassOneError` |
| [`pass_one/tests.rs`](../../../src/ng/run/paralog_filter/pass_one/tests.rs) | 355 lines: 8 tests |

Both subcommands gain `--paralog-fdr` and `--paralog-filter-tag` under their `Advanced` heading,
a refusal variant on their error enum, and the guard that raises it. Seven test-side construction
sites gained the two fields, every one of them named by the compiler.

## Assumptions and recorded deviations

**1. `--paralog-fdr` defaults to `0`, where spec §3.6 says `0.01`.** The spec's default turns the
filter *on*, and the passes that finish it — scoring and writing — are steps C3 and C4. Shipping
§3.6's default now would mean every run parked its records and produced no VCF. **The default
becomes `0.01` at C4**, when the whole path exists; until then the flag is real, its off value is
the one the spec gives it, and nothing is silently different.

**2. A non-zero target is refused rather than run half-way.** Pass one is built and a run could
park every record — and then have nowhere to go. Refusing means the flag's plumbing is exercised
end to end (the CLI parses it, the run reads it, the error names it) while no run can reach a
state where the operator has a spill and no calls. **The refusal happens before the writer is
opened**, so a refused run leaves nothing behind, which is what the test asserts.

**3. The sink is a type, not a closure at each call site.** Two subcommands make the same choice,
and the choice has to be identical in both — `CalledRecordSink` is where "identical" is
enforceable, and C4 changes which variant is built rather than touching either call site again.
Its off variant holds the `VcfWriter` and calls `write_record` with nothing in between.

**The sink closes itself**, rather than handing the writer back to be finished. An earlier version
pattern-matched the off variant out and declared the other `unreachable!` — sound only because a
guard ninety lines earlier had refused it, which is a panic waiting for the step that moves that
guard. `CalledRecordSink::finish` consumes the sink and does the right thing per variant.

**Corrected at review:** the same version also built the writer *before* choosing the sink, so
this module's central promise — with the filter on the VCF is never opened — was a property of the
type that the wiring could not honour. The writer is now built inside the choice.

**4. `entry_for` is a free function, not a method.** It is the whole of the record → entry rule
and it needs no state beyond its four arguments, so it is testable without a sink, a file, or a
run. That is where seven of the eight tests point.

## Tests

**17** — twelve of the entry rule and the sink, and five of the two CLI guards. Nine as first written; the review added eight.

**The one that would have been silently wrong:** `a_left_padded_deletion_is_parked_at_the_position_it_is_written_at`.
A deletion that pads left is *written* one base before its span, and pass three re-runs the
writer's ordering check from the entry's head fields — so an entry carrying the span start would
order the file against a base the file does not contain, and the two would disagree by exactly one
on every padded deletion. The entry's position comes from the same function the `POS` column does,
and the test asserts the two agree by reading `POS` out of the encoded line.

**On the rule from §3.2, at the point it is applied:** a repeat tract is parked with windows and
no read counts; a multiallelic site is parked with its alternatives *summed* (2 reference, 3 and 4
alternative → `ref_reads = 2, alt_reads = 7`); a record with no alternative parks a zero.

**On what must not be scrambled:** the windows reach the entry in the run's sample order, checked
by giving three samples three different depths *and* three different read counts and asserting
both arrive together — a swap that moved only the windows would pass a test that checked only
depths.

**On the head fields as a whole:** for a padded deletion and for a tract, the place pass three
would rebuild from the entry is compared against the position the encoded line carries and against
the record's own tract flag.

**On the sink:** parking a record does not create the VCF, and does create the spill.

**And on the guard:** a non-zero target is refused, with the target named, and no output file
appears. Mutated by making the guard unreachable: the test fails.

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::run::paralog_filter::pass_one"
    → test result: ok. 14 passed; 0 failed; 0 ignored; 6581 filtered out

    cargo test --all-features --lib "call_from_psps::tests::a_paralog_target"
    → test result: ok. 1 passed; 0 failed

    cargo test --all-features --lib --bins --tests
    → 6,582 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 9 errors: 3 needless_lifetimes in cohort_merge, 6 in src/ng/window_coverage/
      (counted by grouping every error kind, not the three the baseline had — see below)

    cargo fmt --check
    → dirty on 9 files, none of them this step's

The integration failure and the three `cohort_merge` lints are `main`'s at `a33ada0f`; the six
`window_coverage` lints and the nine files arrive with the merged branch and fire there
identically. The gate holds, and the lib count moves 6,565 → 6,582.

**As first written this step reported that gate wrongly**, and the way it did is worth more than
the correction: it added two `variable does not need to be mutable` errors on its own lines, and
the command checking the gate grepped for the three lint kinds the baseline already had — so a new
kind was invisible by construction while the count kept matching. Both lints are gone, and the
figure above comes from a command that groups every kind.

## What this step does *not* prove, and why the plan asked for it

The plan asks for byte-identity **on the run fixtures and the tomato slice**, and this step does
not run them. What it offers instead is that the off path's only change is a delegation:
`accept`'s off arm is `writer.write_record(record)` and nothing else, so the bytes cannot differ
unless that one line is wrong — and the writer's own suite already pins what `write_record`
produces. **That is an argument, not a measurement**, and the difference is worth stating: a
fixture run would catch a mistake in the wiring *around* the sink that reading the diff might
miss. It belongs with C5's oracles, which need the on-path to exist before any of them can run.

**What it cannot establish is the on-path**, which has no output to compare until pass three
exists. The four oracles in step C5 are where that lands — off identical, unflagged identical,
drop equals tag minus the tagged, and the two modes identical — and C5 is after C4 for exactly
this reason.

## Tradeoffs and follow-ups

- **Nothing reads the spill yet.** `SpillingSink` parks entries and can hand back the file; the
  scoring pass is C3.
- **The default moves at C4**, with the flag's help text and spec §3.6 agreeing again.
- **`--paralog-filter-tag` is parsed and unused** until C4 gives it a verdict to change. Its help
  text says so rather than leaving an operator to find out.
