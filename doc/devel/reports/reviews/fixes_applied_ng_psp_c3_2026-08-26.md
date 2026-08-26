# Fixes applied — ng psp C3, the body stands on its own

*2026-08-26. Answers [the review](ng_psp_c3_2026-08-26.md) of commit `408f691d`, finding by
finding. Branch `ng-psp-encoding`.*

---

## What changed, in one paragraph

**Both Blockers and all nine Majors are fixed**, and the two Blockers were the same thing: an oracle
that could not fail. The harness now asserts that it built exactly what it was told to build, and
the fixture now carries chain ids, so Milestone E's arrival is covered rather than merely promised.
The sampled property test became an exhaustive one over all 4,096 patterns, which also retired the
fossil seed file the commit had checked in. The record count is 72, up from 70 — fewer tests than
before in one place, because a test that asserted only half its own name is gone.

## The oracle was re-falsified after being changed, three ways

A fix to a test is worth nothing until the test can still fail. All three run in the container, and
each was reverted from a byte-for-byte copy:

| what was injected | before the fixes | after |
|---|---|---|
| a harness that builds nothing | **70 passed, 0 failed** — the Blocker | **67 passed, 5 failed** |
| `as_a_decode_can_return_it` made the identity, i.e. the comparison reaches the chain ids | *(the helper did not exist)* | **66 passed, 6 failed** |
| a body field coded as a step from the previous record | 64 passed, 6 failed | **64 passed, 8 failed** |

The middle one is the check that B2 is really closed: with the stripper removed the comparisons
reach the fixture's chain ids, so the day Milestone E writes them the oracle covers the exception
lists with no test touched.

## Finding by finding

### Blockers

**B1 — the harness never checked what it built.** *Fixed.* `walk_records` (renamed from
`walk_keeping`) asserts `built.is_some() == is_kept(index)` for every record, inside the harness, so
every test that runs through it inherits the check. It also now compares the head the building path
returned against the head the skipping path read, and both byte counts.

**B2 — the fixture carried no chain ids.** *Fixed.* `twelve_records_in_order` puts distinct chain
ids on some observations of the tract records, and every comparison goes through
`as_a_decode_can_return_it`, which strips them because C1's encoder drops them. **When E starts
writing them that helper becomes the identity and can be deleted** — and the measurement above shows
the ids are genuinely in the comparison path today, not merely present in the fixture.

### Majors

**M1 — the fixture's variety was unguarded.** *Fixed.* `the_run_is_the_run_its_doc_describes`
asserts every clause the fixture's doc claims: the record count, spans from one to twelve, **no gap
equal to its own record's span** for all eleven pairs, all three locus kinds, tract payloads that
differ from each other, reference bases that differ between records of the same length, records with
no observations, chain ids present, and at least one record whose head fields each need more than a
byte — that last measured from the encoded heads rather than asserted.

**M2 — the property test sampled about 6 % of the masks.** *Fixed.*
`every_skip_pattern_builds_exactly_the_records_it_keeps` enumerates all 4,096, deterministically.
**That also retires the fossil**: `proptest-regressions/ng/psp/record.txt` is deleted, and this
commit says why — the seed in it was recorded while the *deliberately injected* defect was in the
tree, so it never described a real failure of the shipped code and would have replayed for ever.

**M3 — identical tract payloads.** *Fixed:* a different motif and different flanks on every tract
record, and the fixture test asserts no two are equal.

**M4 — a test asserting half its own name.** *Fixed by deleting it.* Both halves were already
covered: the first is the `"none at all"` pattern of the test above it, and the second is C2's
`a_base_stale_by_one_record_moves_every_record_after_it`, which prices the same fault harder. The
contrast it existed to draw — a body may be skipped, a head may not — is now stated where it belongs,
on `encode_record_body`.

**M5 — reference bases all one byte value.** *Fixed:* they now differ between records of the same
length, and the fixture test requires it.

**M6 — `a_block_of_records` used "block" before blocks exist.** *Fixed:* `twelve_records_in_order`,
with a doc line saying why it is not called a block — a psp block is a span of reference compressed
as one unit, cut on a 100 kb grid, and this is 232 bases in memory. **Third review running to unpick
a reused term**, so the doc says the reason rather than just the name.

**M7 — the record count in three places.** *Fixed:* `RECORDS_IN_THE_RUN`, with the fixture asserting
it and the "only the last" pattern deriving from it.

**M8 — the mask capped the fixture at sixteen records.** *Fixed:* a `u32` mask, and
`const _: () = assert!(RECORDS_IN_THE_RUN < 32);` makes the remaining ceiling a compile error.

**M9 — "the signature is the guarantee".** *Fixed, by making the claim true to what was measured:*
the doc now says the signatures are most of the guarantee and not all of it, that a future edit
could still reach state through a `static` or a thread-local without changing either — which was
measured to compile — and that the tests are what hold the line.

### Minor and Nits — what was taken

Applied: the "every read agreed" arm's observation now really matches the reference it is given, so
the comment and the record agree; the arm that inherits the rich fixture's holed witness asserts its
span is wide enough to hold it, so the fixture stops containing a record no producer could emit; the
300-read builder is extracted and shared with C2's depth test, so the depth this caller is committed
to is written once and the copy that carried no assertions is gone; the compare loop is one helper
used by every test; the harness takes the walk's name, so a refusal says which of six patterns was
running; the scribble runs **all six patterns and two kinds of damage** — `0xff`, which no body can
end on, and a second filling that **is** a decodable body, which catches the one shape the first
cannot: a reader that touches a skipped body, decodes it, and uses what it found; a new test walks
the run under a layout declaring two unknown trailing fields, the configuration the mean-coverage
trap would arrive in; another skips by what the head says rather than by an ordinal, which is the
skip the cohort's first pass actually makes; the bare participles are named; and the doc heading is
a noun phrase again with the one-off tally moved to this report, where its setup is.

**Three claims corrected**, all mine and all about my own work: "a byte the encoder never writes" —
it writes `0xff` eighteen times in this very run, and the true and stronger statement is the one the
inline comment already made, that no body can *end* on it; "one at the top of the depth range" — it
is three; and "it is the only thing that would fail" — six tests do.

**Not taken:** splitting the test module. Three agents' worth of evidence says no: `record.rs` is the
sixth largest of about sixty ng files and all five larger ones keep their tests in-file, 123 ng files
use an in-file `mod tests` against one that does not, and C3 adds no third concern beyond the two the
architecture document assigns.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::record` | **72 passed** |
| `cargo test --lib ng::psp` | **116 passed** |
| `cargo test --lib --bins --tests --examples` | library **4,653 passed**, 14 ignored; every other target green except `examples/ng_generic_loci_dump` (11 failures, all `NotFound { ".../ref.fa.repeats.parquet" }`, pre-existing) |
