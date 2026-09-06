# The hidden-duplication filter — B2: the spill file, and why it cannot outlive its run

**Date:** 2026-09-06
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone B, step B2
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.4, §5
**Branch:** `ng-paralog-filter`

## The answer

**The spill is now a file — `<output>.paralog-spill.tmp`, beside the output — and it goes away
when the run ends, whatever way the run ends.** Success, an error travelling up through `?`, and
a panic unwinding through the caller all remove it, because the removal is a `Drop` rather than a
call: the value going out of scope *is* the cleanup.

That distinction is the step's whole content, and it is measurable. A guard written as a call on
the tidy path — removing the file only where the writer was properly finished — passes the
normal-exit test and **fails the other three**. The spill holds every record's line uncompressed,
so against a `.vcf.gz` output it is several times the output's size — and the runs where leaving
one behind costs the most are exactly the runs a call-shaped guard would miss.

**The file appears on the first record and not before.** A run with the filter off names nothing
and leaves nothing; a run that called no records at all still ends pass one with an empty spill,
which is a different thing from a spill that was never opened, and which lets passes two and
three read an empty stream rather than meet a missing file.

## What was added

| file | what it is |
|---|---|
| [`src/ng/run/paralog_filter/spill_file.rs`](../../../src/ng/run/paralog_filter/spill_file.rs) | 406 lines: `SpillFile`, its three-stage life, `SpillFileError`, and the `Drop` that unlinks |
| [`src/ng/run/paralog_filter/spill_file/tests.rs`](../../../src/ng/run/paralog_filter/spill_file/tests.rs) | 530 lines: 19 tests |

[`spill.rs`](../../../src/ng/run/paralog_filter/spill.rs) gains 44 net lines — the completeness
check below and the sample-count ceiling — and
[`mod.rs`](../../../src/ng/run/paralog_filter/mod.rs) its declaration and re-exports.

**The shape.** `SpillFile::beside(output)` names the file and creates nothing. `append` creates it
on the first entry — with `create_new`, so anything already at the path is named rather than
destroyed — and refuses a record once pass one has ended. `finish_writing` ends pass one: it
flushes, and on a run that appended nothing it creates the file empty. `read` hands back a cursor
from the beginning and may be called more than once, because pass two scores the entries and pass
three writes them. Dropping the value closes the handle and unlinks the file, warning on stderr if
the removal fails.

**Which of the three stages the file is in is a value, not an inference.** It was read off
`writer.is_some()`, and the review found what that costs: a record appended after pass one had
ended reopened the file with `truncate(true)`, so a 228-byte three-record spill became 77 bytes
and one record while `entries_written` said 4 — and the read then reported a *lost tail* on a file
that had lost its head.

## Assumptions and recorded deviations

**0. The review found six Major defects and this section is the state after them.** The largest
was that an append after pass one had ended silently emptied the whole spill; the file's stage is
now a value rather than an inference, and creation refuses to overwrite anything already at the
path. The
[review](../reviews/ng_paralog_filter_b2_2026-09-06.md) and the
[fixes applied](../reviews/fixes_applied_ng_paralog_filter_b2_2026-09-06.md) carry the rest,
including three numbers of mine they corrected.

**1. The completeness check the B1 review deferred to this step is built here.** B1's review found
that a spill losing its tail on an entry boundary reads back as a complete, shorter spill — the
bytes cannot tell the two apart, because a truncated file *is* the prefix of a longer one — and
that pass three would then write a VCF short by its last records with no panic and no message.
The count exists on the writer and nothing carried it across. It does now: `SpillReader::new` **takes the count as an
argument** — it was a separate `expecting` call, which the review pointed out is one forgotten
line away from the failure this closes — and a file that ends holding a different number gives
`SpillError::TheWrongNumberOfRecords { expected, read }`, in either direction. **This is more than the plan's
step B2 asks for**, and it is here because this is the first step where the number and the reader
are in the same place.

**2. `read` refuses a spill that pass one has not finished**, rather than returning what happens
to have reached the disk. At 64 KiB of buffering that is usually nothing at all, so the failure
would be a spill silently short by everything. `SpillFileError::PassOneHasNotEnded` names it as a
caller's ordering mistake, and the same variant covers reading a spill that was never opened —
the answer there is "pass one has not finished", not "the file is missing".

**3. `finish_writing` creates the file on a run that appended nothing.** The alternative is for
`read` to hand back an empty stream without a file behind it, and there is no portable empty
`File` to build that from. Creating it at the end of pass one is both simpler and truer: pass one
ran, so there is a spill, and it holds no records. A run with the filter off never calls
`finish_writing` and leaves nothing on disk.

**4. A failed removal warns on stderr and cannot do more.** A `Drop` running while a panic unwinds
cannot return an error and must not panic itself, so the warning is skipped while
`std::thread::panicking()` — a second message there buries the first. The crate already answers
this the same way in [`reference_info`](../../../src/ng/reference_info.rs)'s `VerificationHandle`.
Silence was the first draft's choice and the review was right that it is wrong: what is leaked is
a file several times the size of a compressed output.

**5. The path is appended to the whole output path, not substituted into it.** `cohort.vcf`
becomes `cohort.vcf.paralog-spill.tmp`, which is the convention
[`vcf::writer`](../../../src/ng/vcf/writer.rs) already uses for `<output>.tmp`, and which gives one
predictable answer for an output with no extension, several, or a name that looks like a
directory.

**6. Both directions are buffered at 64 KiB**, matching the VCF writer. B1's review noted that
`append` issues one `write_all` per entry, so an unbuffered sink would cost a write syscall per
called record; the reader's `BufRead` bound says the same for its byte-at-a-time varints. This is
the step that chooses the sink, so it is the step that pays the requirement.

## Tests

**19 under `ng::run::paralog_filter::spill_file::`**, seven of them from this step's review; see
the [fixes applied](../reviews/fixes_applied_ng_paralog_filter_b2_2026-09-06.md).

**Four are the exits**, and they are separate tests because they are three different mechanisms
plus one variation: a value going out of scope at the end of a block; a value going out of scope
because `?` returned a `RunError` early; a value going out of scope because a panic is unwinding
past it; and a run that ends with the writer still open, never having flushed.

The other fifteen: the path is the output's with a suffix added, over three shapes of output
name; the file does not exist until the first record and does after it; a spill named and never
written to creates nothing **and removes nothing**, including a file it did not create; something
already at the path is refused rather than overwritten, and so is a second spill for one output;
an append after pass one has ended is refused and the file keeps its bytes; a record the codec
refused is not counted; the records come back from a real file bit for bit, with the absent sample
still absent; the file can be read twice, the second cursor starting over; a run that called
nothing reads an empty spill; finishing pass one twice leaves the records where they are; reading
before pass one has ended is refused, both before the file is opened and before it is flushed; a
spill whose tail was lost is refused rather than read short; a read failure can be given the file
it happened on; and a spill that cannot be created names the path and what the filesystem said.

**Three mutations, each killed, each reverted from a byte-compared backup.** Run as
`cargo test --all-features --lib "ng::run::paralog_filter::"` in the container:

| mutation | what happened |
|---|---|
| the `Drop` deleted outright | 4 of 19 fail — all four exits |
| the guard kept but fired only where pass one was properly finished — *a guard on the tidy path* | **3** of 19 fail, and **`the_file_is_gone_when_the_run_ends_normally` passes**: the shape that looks correct on the happy path and leaves the file on every other |
| an append after pass one reopening the file with `truncate`, as before the review's fix | 3 of 51 fail |

**The second mutation found a defect in the tests as well as proving them.** With the file left
behind by one failing test, `the_file_does_not_exist_until_the_first_record_is_written` failed on
the *next* run rather than on its own merits — a test asserting a file is absent was reading what
an earlier run had left. `an_output_path` now removes any stale spill before handing the path
over. **The count reported here was wrong because of it**: the first draft of this report said the
mutation fails four tests, and four was what that contaminated run showed. Re-measured after the
fix, it is **three** — which is the number that carries the point, since the fourth was never one
of the exits.

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::run::paralog_filter::"
    → test result: ok. 51 passed; 0 failed; 0 ignored; 6373 filtered out

    cargo test --all-features --lib --bins --tests
    → 6409 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/run/paralog_filter/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f`. The gate —
no failure `main` does not already have, and the `fmt` and `clippy` failure sets no larger than
`main`'s — holds, and the lib count moves 6,387 → 6,409.

**`cargo fmt` reformats `main`'s nine files as a side effect and they were reverted with
`git checkout --` before staging.**

## Tradeoffs and follow-ups

- **The scratch directories the tests create under `tmp/paralog_spill_file_tests/` are left in
  place**, empty, because each test removes its own file and nothing removes the directories. They
  are inside the project's ignored `tmp/`, per `CLAUDE.md`.
- **How much larger the spill is than a compressed output is not measured**, only its direction.
  Step D2 is the first run that could measure it.
- **The spill's disk cost is reported nowhere.** Spec §3.5's run-report list has no line for it,
  so adding one is a spec question; raised at Checkpoint B.
- **A spill left by a killed run now blocks the next run of the same output** until it is removed,
  because creation refuses to overwrite it. That is the right way round — the alternative destroys
  a finished run's evidence or a live run's records — but it is a behaviour an operator meets.
- **A run killed from outside still leaves the file.** Nothing inside a process runs when the
  process is killed; the spec says so (§3.4) and `<output>.tmp` behaves the same way.
- **Nothing appends to a `SpillFile` yet.** The sink that fills it is step C2, which is also where
  a `SpillFileError` becomes the `RunError` naming the file that spec §5 asks for.
