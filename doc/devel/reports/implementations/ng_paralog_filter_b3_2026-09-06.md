# The hidden-duplication filter — B3: writing a line whose record is gone, and moving two columns of it

**Date:** 2026-09-06
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone B, step B3 — *own commit, do not bundle*
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.5, §6 traps 5–7
**Branch:** `ng-paralog-filter`

## The answer

**Pass three can now write a record whose record no longer exists, and change two columns of it
without touching the rest.**

Two pieces, and each closes one of the ways this could go quietly wrong.

**The writer gains a way in for a line, and it is the same way in the record takes.**
`write_record` is now `write_line` with the line built first, so there is one ordering check
rather than two that could drift. That matters because the filter parks finished lines on a spill
and writes them back a pass later, by which time the record is gone — and a VCF whose records run
backwards is not a VCF, which no consumer would notice until it indexed one (spec §6 trap 6).
Measured: deleting the check from `write_line` fails **nine** of the writer's twenty tests, and
**seven** of the nine pre-date this step — six calling `write_record` directly and one calling
`write_stream`. There is nowhere left to bypass.

**The patch keeps the line's bytes and splices two columns into them.** It does not rebuild the
line from its parts, which is what makes the filter's standing oracle checkable rather than hoped
for: a run with the filter on at a target no record reaches, with the two `INFO` keys stripped,
must equal the run with the filter off **exactly** (spec §10). A rebuild would have to reproduce
every column's spelling; this keeps them.

## What was added

| file | what it is |
|---|---|
| [`src/ng/run/paralog_filter/patch.rs`](../../../src/ng/run/paralog_filter/patch.rs) | 187 lines: `rewrite_filter_and_info` and `LinePatchError` |
| [`src/ng/run/paralog_filter/patch/tests.rs`](../../../src/ng/run/paralog_filter/patch/tests.rs) | 741 lines: 23 tests, of which 4 encode real records and 2 are properties |

[`vcf/writer.rs`](../../../src/ng/vcf/writer.rs) gains 57 net lines — `write_line`, `RecordPlace`
made public and carrying a `GenomePosition`, and `place_of` reading the padding rule from
`encode::written_position` rather than a second copy of it — with 195 net lines of tests beside it.
[`paralog_filter/mod.rs`](../../../src/ng/run/paralog_filter/mod.rs) gains the declaration and
`From<&SpillEntry> for RecordPlace`; [`vcf/mod.rs`](../../../src/ng/vcf/mod.rs) re-exports
`RecordPlace`; [`vcf/encode.rs`](../../../src/ng/vcf/encode.rs) makes `written_position`
`pub(crate)`.

**The line and test counts above are the step as it stands after its review**, whose fixes are
folded into this commit; as first written `patch.rs` was 139 lines with 13 tests and `writer.rs`
gained 28.

## Assumptions and recorded deviations

**1. `write_record` is expressed in terms of `write_line`, which the plan does not ask for.** The
plan asks for the entry point to run "the same ordering check `write_record` runs". Two call sites
running the same check is a thing that can drift; one check that both paths go through cannot.
The cost is that `write_record` no longer holds the write itself, and the benefit is what the
mutation above measures — the check cannot be removed from the line path alone.

**2. `RecordPlace` becomes public.** It carries what a caller holding a line and no record has to
hand over. The alternative — three parameters — spells the same thing without a name for it, and
`place_of` would still have to return something. **Corrected at review:** as first written it held
a `ContigId` beside a bare `u64`, which is *not* the shape `SpillEntry` carries (`position:
Position`), so pass three would have had to unwrap at the one call site where a place that
disagrees with its line is accepted. It now carries a `GenomePosition`, and
`From<&SpillEntry> for RecordPlace` is how pass three builds one.

**3. The patch lives in `paralog_filter/`, not in `vcf/`.** Rewriting a line's seventh and eighth
columns is VCF arithmetic and could sit beside the encoder, but its rules are the filter's: which
`FILTER` values may be joined to, and which `INFO` values are replaced, come from spec §3.5. Put
in `vcf/` it would be surface with one caller in another module.

**4. The "nothing to add" case goes through the splice rather than around it.** A short circuit
returning the input would make the plan's first test — a patched line with nothing to add is
byte-identical to its input — pass without saying anything about the splice. Going through it
means that test *is* the proof the rebuild is faithful.

**5. A `FILTER` of `.` is replaced, like `PASS`.** Spec §3.5 names the join rule for a record that
"already carried a filter" and ng's own encoder never writes `.` there — `FilterVerdict` has five
values and all are words. But `.` is VCF's *no filters applied*, so joining to it would produce
`.;hiddenParalog`, which says two contradictory things. Both `.` and an empty column are treated
as nothing to join to.

**6. An empty `INFO` column is replaced, like `.`.** The spec names `.`; the empty column is
handled the same way, since appending to either would put a leading `;` on it.

**Corrected at review:** this was justified by saying ng's encoder writes an empty `INFO` for a
record with no annotations, "because `info_column` joins an empty field list". **It does not and
cannot** — `encode.rs:218-219` pushes `AN=` and `DP=` before any conditional field, so the list is
never empty. Measured over every shape ng emits, the thinnest `INFO` is `AN=0;DP=0`. Both spellings
are handled defensively, for a line from somewhere other than ng's own encoder; neither is
reachable, and the tests that pin them now say so in their names.

**7. A line that does not split into nine pieces is refused, not patched.** The nine are the eight
fixed columns, then `FORMAT` and the sample columns together. A shorter line did not come from the
encoder, and rewriting its seventh piece would be rewriting something other than `FILTER`.

**Corrected at review:** as first written, the variant's doc said "a record's line always has at
least ten", the `# Errors` section said "the eight columns a VCF record has before `FORMAT`", and
the code refused at fewer than nine — three rules, no two alike, and a nine-column line was
accepted while the doc said it was not. The invariant the doc stated is true —
`VcfRecord::new` refuses a record with no sample columns, so the encoder cannot produce nine — but
it is not what this guard is for. All four statements now say nine, and the doc says why the tenth
is not checked.

## Tests

**23 in `patch`, 7 new in `vcf::writer`** (which now has 22). Thirteen and five of those are the
step as first written; the review added ten and two.

**The patch, on what must not move:** a line with nothing to add comes back byte for byte; with
both a filter and two `INFO` fields added, **every column but the seventh and eighth is compared
and must be unchanged**; the filter and the `INFO` are independent, so a filter-only patch leaves
the `INFO` alone and the reverse; and a record with one sample and one with **a thousand** patch
the same two columns, with all thousand sample columns compared. A run of tabs inside the sample
columns survives, because the split stops after the eighth.

**On the two rules that are easy to reverse:** a passing record carries the filter in place of
`PASS`; a record already on `EMNoConv` carries `EMNoConv;hiddenParalog`; a `FILTER` of `.` or empty
is replaced rather than joined to; an `INFO` of `.` and an empty `INFO` are both replaced rather
than appended to; and the fields are appended in the order they are given. **The three tests of
unreachable spellings say so in their names** (`…_though_ng_never_writes_one`), and the reachable
half of trap 7 — a real verdict joined to — is pinned in `encoded_lines` on lines the encoder made.

**On lines the encoder actually made** (`encoded_lines`, added at review): twelve shapes — every
`FilterVerdict` spelling, both `FORMAT` strings, both padding sides, `ALT .`, a no-called sample, a
tract whose `REPCN` is `.`, a record on the second contig, and a three-thousand-sample cohort —
are built as `VcfRecord`s, encoded through `record_line`, and required to come back byte for byte
with nothing to add. A second test says only columns seven and eight move; a third says what those
two *become*, which is the half that no test had. **This is spec §10's oracle at record scale**,
and the committed hand-written fixtures were not the encoder's shape: `AF` at two decimals where
the encoder fixes six, and a `FORMAT` of `GT:GQ:AD` where it writes `GT:GQ:DP:AD`.

**Two properties** over 10 to 40 arbitrary columns: nothing to add is the identity, and with
something to add only the seventh and eighth move.

**On what is refused:** a seven-column line, an eight-column line, and an empty one, each naming
how many columns it found.

**The writer:** a line written without its record reaches the file exactly; **writing a record and
writing its line produce the same bytes**; a line that runs backwards is refused as a record would
be; a line may take the one legal tie and no other; and a line counts towards `records_written`.

**Three mutations, each killed, each reverted from a byte-compared backup.** Run in the container:

| mutation | what happened |
|---|---|
| the ordering check deleted from `write_line` | **9 of 20** writer tests fail — **seven of them written before this step**, because both paths now go through the one check |
| a pre-existing `FILTER` replaced rather than joined — spec §6 trap 7 | 1 of 13 fails: `a_record_that_already_carries_a_filter_keeps_it_and_gains_the_new_one` |
| a `.` `INFO` appended to rather than replaced — spec §3.5 | 1 of 13 fails: `an_info_column_that_says_nothing_is_replaced_rather_than_appended_to` (renamed at review) |

**Three more mutations survived the thirteen and were closed at review**: dropping the empty-`FILTER`
guard, deleting `append_info_column`'s early return (which turned a `.` `INFO` with nothing to add
into an empty column, breaking the identity this module promises), and advancing the writer's
position and record count before the ordering check. Each now has a test.

**The first is the step's structural claim, measured.** Nine failures where a separate check would
have given one is the difference between "the line path also checks" and "there is one check".

## Validation

Run in the dev container, on this tree, **after the review's fixes were folded in**.

    cargo test --all-features --lib "ng::run::paralog_filter::patch"
    → test result: ok. 23 passed; 0 failed; 0 ignored; 6431 filtered out

    cargo test --all-features --lib "ng::vcf::writer"
    → test result: ok. 22 passed; 0 failed; 0 ignored; 6432 filtered out

    cargo test --all-features --lib --bins --tests
    → 6439 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them this step's

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f`. The gate —
no failure `main` does not already have, and the `fmt` and `clippy` failure sets no larger than
`main`'s — holds, and the lib count moves 6,409 → 6,439 (6,427 as first written, plus the review's
twelve).

## Tradeoffs and follow-ups

- **Nothing calls the patch yet.** Pass three is step C4, which is also where the filter id and
  the two `INFO` field names come from — this step takes them as arguments and knows neither.
- **The patch allocates a fresh `Vec` per line.** At a cohort's line length that is one allocation
  a record in pass three — about 33 kB at three thousand samples — and whether it is worth a reused
  buffer is what step D2's pass shares say. **The owned return means adopting one later changes
  both this signature and C4's call site**, which the review raised and this step deferred: a
  `…_into(&mut Vec<u8>, …)` form would let D2's answer land without touching a caller, but it adds
  a second entry point with no caller today.
- **`write_line` does not check the line's shape**, deliberately: it owns the order and the bytes'
  destination, and the caller owns what the line says. **The patch refuses one malformedness of
  three:** too few columns. A `\n` in the line, or a `\t` in an added value, passes through and
  silently reshapes the file — stated on `rewrite_filter_and_info` rather than checked, because
  neither is reachable from `record_line` and pass three's added values are C4 constants.
- **Where the module lives is open.** Deviation 3 argued it belongs in `paralog_filter/` because
  its rules are the filter's. The review inventoried them and all five are VCF grammar — the tab
  separator, the column positions, `.`, `PASS`, `;` — with the filter id and the two keys arriving
  as arguments. Two of the five are now taken from `src/ng/vcf` rather than re-spelled; whether the
  file itself moves is a Checkpoint B question, since it reverses a recorded deviation.
- **The error names no record**, which is untraceable on a file of millions. `rewrite_filter_and_info`
  is handed bytes and nothing else, so wrapping it with the entry's contig and position is C4's
  job; the doc now says so.
