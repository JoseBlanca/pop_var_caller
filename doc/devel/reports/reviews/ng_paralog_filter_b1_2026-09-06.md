# Code Review: ng_paralog_filter_b1

**Date:** 2026-09-06
**Reviewer:** rust-code-review skill (orchestrator), eight sub-agents in isolated worktrees
**Scope:** step B1 of the hidden-duplication filter plan — the spill entry and its codec
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the working-tree diff of commit `2ad2a92b` on branch
  `ng-paralog-filter` — three new files plus one added `pub mod` line.
- **Reviewed against:** `2ad2a92b`. Every sub-agent detached its own worktree to that commit and
  confirmed two branch-only files were present before starting.
- **In-scope files:**
  - [spill.rs](../../../../src/ng/run/paralog_filter/spill.rs)
  - [tests.rs](../../../../src/ng/run/paralog_filter/spill/tests.rs)
  - [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs)
  - [mod.rs](../../../../src/ng/run/mod.rs) — the one added line
  - [ng_paralog_filter_b1_2026-09-06.md](../implementations/ng_paralog_filter_b1_2026-09-06.md) — its quantitative claims
- **Deliberately out of scope:** `src/ng/paralog/` (Milestone A, reviewed and fixed);
  `src/paralog/`, `src/var_calling/`, `src/sample_summary/` (production, frozen);
  `src/psp/varint.rs` (called as-is).
- **Categories dispatched:** reliability, errors, extras (the code is a parser producing stable
  output), refactor_safety (a type is a stand-in for another branch's), naming, idiomatic,
  smells + module_structure (a new module), defaults + the diff's own numbers. `tooling` was
  skipped — `Cargo.toml` is untouched; `unsafe_concurrency` was skipped — no `unsafe`, no
  threads, no shared state.

### 2. Verdict

**Request-changes.** One Blocker: the comparator that nine round-trip tests rest on cannot be
shown to fail, and hard-wiring it to `true` leaves all seventeen tests green.

### 3. Execution status

| command | result |
|---|---|
| `cargo test --all-features --lib "ng::run::paralog_filter::"` | `ok. 17 passed; 0 failed; 0 ignored; 6373 filtered out` |
| `cargo test --all-features --lib --bins --tests` | 6,375 lib tests passed; one integration test failed, pre-existing on `main` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 3 errors, all `needless_lifetimes` in `cohort_merge/`, pre-existing |
| `cargo fmt --check` | dirty on 9 files, none of them this commit's, pre-existing |

`--all-targets` was not run: `examples/ng_candidate_selection_probe.rs` has not compiled on
`main` since `a33ada0f`, so the command measures nothing. `cargo doc --no-deps` and
`cargo audit` were not run — the first is red on `main` for 11 unresolved intra-doc links in
other modules, the second has no offline database in the container.

**Findings labelled "Needs verification": 0.** Every finding below was produced by mutation or
by a probe run in the reviewer's own worktree, and the outputs are quoted in the per-category
files under [tmp/review_2026-09-06_ng_paralog_filter_b1/](../../../../tmp/review_2026-09-06_ng_paralog_filter_b1/).

### 4. Open questions and assumptions

1. **Does `SpilledSample` keep the spec's flat shape?** Spec §3.7 gives `SpillEntry` two
   independent `bool`s and every `SpilledSample` an `alt_reads`, while §3.2 describes one
   either/or: a biallelic SNP is scored on both allele counts, every other record — repeat
   tracts included — on coverage alone. Two agents independently proposed a sum type. That
   changes the spec's type sketch, so it is the owner's call. Affects **M10**.
2. **Who owns the spill's completeness?** Nothing in the file says it is finished, so a spill
   that loses its tail on an entry boundary reads back as a shorter valid spill. The count
   exists on the writer (`entries_written`) and nothing carries it to the reader. B2 owns the
   file's lifecycle and is the natural place to thread it. Affects **M9**.
3. **How far does the error context go before the run wiring adds its own?** Spec §5 says a
   spill failure becomes a `RunError` naming the file, which is C2–C4's. The record ordinal is
   also the caller's; the sample index is only this module's. Affects **M6**, **M8**.

### 5. Top 3 priorities

1. **B1** — the round-trip comparator cannot fail, so nine tests assert only that decoding did
   not panic.
2. **M1** — a corrupt line length pulls the rest of the file into one `Vec`, measured at 16 MB.
3. **M4** — no fixture has a length prefix at or above 128, which is the regime every real
   record is in.

### 6. Findings

#### Blocker

- **B1: [tests.rs:75](../../../../src/ng/run/paralog_filter/spill/tests.rs#L75) — `same_bits` is what nine round-trip tests assert, and nothing shows it can return `false`**
- **Categories:** reliability
- **Confidence:** High — measured.
- **Problem:** hard-wiring `same_bits` to `return true` leaves **all 17 tests green**. Four of
  the nine round-trip tests then assert nothing but *the decode did not panic* —
  `a_record_whose_info_column_is_a_dot_round_trips`, `an_empty_line_round_trips`,
  `values_at_the_top_of_their_fields_round_trip` and
  `a_stream_of_records_comes_back_in_the_order_it_was_written` have no assertion outside
  `assert_round_trips`. The size of the loss is measurable: the author's own headline mutation
  — a decoded `NaN` turned into `0.0` — kills 7 of 17 with the comparator intact and **2 of 17**
  with it neutered, so five of those seven kills are the comparator's, and the comparator itself
  is never exercised on a pair that differs.
- **Why it matters:** `values_at_the_top_of_their_fields_round_trip` is the only test covering
  `ContigId(u32::MAX)`, `Position(u64::MAX)` and `u32::MAX` read counts. A decoder that returned
  `Position(position + 1)`, dropped a large position's high varint bytes, or swapped `ref_reads`
  with `alt_reads` produces wrong results without panicking, and the comparator is the only thing
  between that and a green suite.
- **Fix:** a negative test of the comparator — one pair per field, including **an absent window
  against a zeroed one**, which is spec §6 trap 4 stated as an assertion about the comparator
  rather than about the codec. Full body in
  [reliability.md](../../../../tmp/review_2026-09-06_ng_paralog_filter_b1/reliability.md).

#### Major

- **M1: [spill.rs:311-319](../../../../src/ng/run/paralog_filter/spill.rs#L311-L319) — a corrupt line length pulls the rest of the file into one `Vec`**
- **Categories:** extras
- **Confidence:** High — measured.
- **Problem:** `line_length` is decoded with no ceiling and handed to
  `take(line_length).read_to_end(&mut line)`. `read_to_end` appends every byte available before
  the `read != line_length` check runs. Fed a header claiming `u64::MAX` followed by a 16 MiB
  payload, the decoder pulled **16,777,230 of 16,777,216 payload bytes into one `Vec`** before
  erroring. Two other agents probed the same field with a file holding *no* payload and saw no
  spike — correctly, and that is the whole mechanism: the ceiling is the bytes that follow, and
  on a real spill those are the rest of the file.
- **Why it matters:** spec §5 requires one entry in hand at a time, and §3.4 sizes the spill at
  about the VCF's size — gigabytes on a whole genome. One corrupt length byte turns a decode
  error into an out-of-memory kill of a run that has already finished calling. The neighbouring
  field is bounded and its doc says why; this one is not.
- **Fix:** refuse a length no VCF line can reach before reading it. A line at 3,000 samples is a
  few hundred kilobytes, so a megabyte ceiling is generous. Add a test asserting the variant.

- **M2: [tests.rs:75](../../../../src/ng/run/paralog_filter/spill/tests.rs#L75) — `same_bits` reads its fields by access, so a field either struct gains drops out of every round-trip test**
- **Categories:** refactor_safety
- **Confidence:** High — demonstrated twice.
- **Problem:** the codec itself is safe: adding a third field to `WindowCoverage` gives five
  compile errors, `#[non_exhaustive]` on top of it gives the same, and a rename gives `E0026`,
  `E0027` and `E0560`. But the oracle is not destructured. Adding `scoring_weight: u32` to
  `SpilledSample`, encoding it correctly and decoding it as `0` where every sample was written
  as `7` leaves **17 passed, 0 failed**; the same on `SpillEntry` with a `bool` written `true`
  and decoded `false` leaves **17 passed, 0 failed**. The compiler forces the decoder to
  *mention* a new field; nothing forces it to decode it *correctly*.
- **Why it matters:** the pinned-bytes test fails in both experiments, but on the *encoder's*
  bytes — the half the destructure already protects — and the coder's correct response is to
  update the pinned list, after which it is green with the decoder still wrong. The four numbers
  per sample are the scorer's whole input.
- **Fix:** destructure both sides with no `..`, exactly as `LocusWindowCoverage`'s `PartialEq`
  and `cohort_merge`'s `render` do — the two precedents this module's own header cites. Re-adding
  the field then fails to compile inside the oracle.

- **M3: [spill.rs:244-262](../../../../src/ng/run/paralog_filter/spill.rs#L244-L262) — a decode error is not terminal: the reader resynchronises and hands back a plausible entry**
- **Categories:** reliability, errors, idiomatic — convergent
- **Confidence:** High — measured.
- **Problem:** `next_entry` returns `Some(Err(..))` and leaves the stream wherever the failure
  stopped it; the next call decodes from that offset. On a file of three bytes that fail on
  their flag byte followed by one whole valid entry, the second call returned
  `Some(Ok(SpillEntry { contig: ContigId(1), position: Position(300), … }))` — a well-formed
  record decoded from the middle of a corrupt file. Separately, a source whose `fill_buf` keeps
  failing yielded 1,000 `Some(Err(…))` in a row and never `None`. The existing offset-walk test
  calls `next_entry()` exactly once per truncation, so it cannot see either.
- **Why it matters:** the module header makes the complete-versus-truncated distinction the
  reason the reader exists, and the `Iterator` impl undoes it. At C2 a consumer that logs and
  continues writes VCF records decoded from a resynchronised offset, or hangs. Three sibling
  fallible iterators in this crate latch — `tandem_repeat.rs:986`, `run/walker.rs:207`,
  `locus_generation/mod.rs:1142`.
- **Fix:** a `failed` latch set by the first failure, plus `impl FusedIterator`, plus the test.

- **M4: [tests.rs](../../../../src/ng/run/paralog_filter/spill/tests.rs) — no fixture has a line or a sample count needing more than one varint byte, which is every real record**
- **Categories:** reliability
- **Confidence:** High — measured.
- **Problem:** reading the line length as a **single byte** instead of a varint passes all 17
  tests. Every fixture line is short — the smallest is 2 bytes, the longest 78 — and every
  `per_sample` is 0, 1 or 2 entries, so both length prefixes always fit one byte. The same
  blindness covers the sample count.
- **Why it matters:** this is the regime the codec will be in on every run and never in a test.
  A 63-sample tomato line runs to several hundred bytes; at the thousand-sample end of the range
  `design_principles.md` §0 commits to, a line is kilobytes and the sample count needs two varint
  bytes. The one case the fixtures cover — a length below 128 — is the one that will never occur.
- **Fix:** two round trips in the regime the caller runs in: a 300-byte line, and a
  thousand-sample cohort. Both verified correct on the reviewed tree, so the gap is coverage,
  not behaviour.

- **M5: [spill.rs:121](../../../../src/ng/run/paralog_filter/spill.rs#L121) — "every variant names the field it failed on" has a regression test for 4 of its 11 labels**
- **Categories:** reliability, refactor_safety — convergent
- **Confidence:** High — measured.
- **Problem:** five separate label mutations survive all 17 tests: `"position"` → `"contig"`,
  `"is_biallelic_snp"` → `"is_repeat_tract"`, `"sample count"` → `"line length"`,
  `"mean_depth"` → `"gc_fraction"`, `"alt_reads"` → `"ref_reads"`. Only `"contig"`,
  `"is_repeat_tract"`, `"line"` and `"gc_fraction"` are asserted anywhere. The test that walks
  every truncation offset asserts `matches!(outcome, Some(Err(_)))` and nothing more — ask which
  wrong implementations satisfy it and the answer is *every implementation that names every field
  `"contig"`*. Separately, the labels are `&'static str` literals: renaming `mean_depth` gave
  four compile errors and **none was the `"mean_depth"` literal**.
- **Why it matters:** the module states its own reason — a spill is written and read inside one
  process, so when it fails the field name is the entire diagnostic, and seven of eleven are
  unpinned.
- **Fix:** strengthen the offset walk to assert the expected label at every cut, which covers all
  eleven in one test; add the missing `is_biallelic_snp` flag-byte test; and hold the labels in
  one place so a rename has a single site.

- **M6: [spill.rs:129](../../../../src/ng/run/paralog_filter/spill.rs#L129) — `Io` is one mechanism-named variant standing in for six operations, and cannot tell a failed write from a failed read**
- **Categories:** errors
- **Confidence:** High — probed.
- **Problem:** `#[from] io::Error` is generated for `write_all`, `flush`, the boundary
  `fill_buf`, `read_to_end`, `read_f32_bits`'s non-EOF arm and `read_byte`, so provenance is
  gone; the message says "could not be read **or** written". A pass-one write failure and a
  pass-two read failure render as the same line. This contradicts the type's own doc four lines
  above it, and the house style next door: `VcfWriteError` splits the identical `io::Error` into
  five operation-named variants.
- **Fix:** drop `#[from]`, add operation-named variants, and give the read arm the field.

- **M7: [spill.rs:130](../../../../src/ng/run/paralog_filter/spill.rs#L130) — `SpillError::Io` has no test, and a `finish` that never flushes passes all 17**
- **Categories:** reliability
- **Confidence:** High — measured.
- **Problem:** deleting `self.sink.flush()?` from `finish` survives all 17 tests, because every
  fixture sink is a `Vec<u8>` whose flush cannot fail. Nothing constructs a stream that refuses a
  write or a flush, so the variant a real run on a full disk hits is never reached.
- **Why it matters:** a `finish` that swallowed the flush error would report a complete spill for
  a truncated one, and pass two would score a cohort missing its last records — a wrong run-wide
  cut, no panic, no message.
- **Fix:** a sink that refuses everything, and two tests over it.

- **M8: [spill.rs:323](../../../../src/ng/run/paralog_filter/spill.rs#L323) — a failure inside the per-sample array names the field but not the sample**
- **Categories:** errors
- **Confidence:** High — probed.
- **Problem:** `for _ in 0..sample_count` discards the index. Cutting a three-sample entry inside
  sample 1 and inside sample 3 gives byte-identical messages. At spec §4's three-thousand-sample
  cohort the message points at three thousand places at once.
- **Fix:** carry the index. The record ordinal is the caller's to add (C2–C4); the sample index
  is only this module's.

- **M9: [spill.rs:209-213](../../../../src/ng/run/paralog_filter/spill.rs#L209-L213) — nothing says the file is complete, so a spill that loses its tail reads back as a shorter valid spill**
- **Categories:** extras, errors, naming — convergent
- **Confidence:** High on the behaviour (measured); Medium on the impact, because no caller
  exists yet.
- **Problem:** `finish`'s doc says consuming the writer means "a forgotten flush shows up as a
  short file at the reader rather than as an entry that quietly never arrived". It does not.
  `SpillWriter` has no `Drop`; taking `self` by value prevents *use after finish* and nothing
  else, and `Drop for BufWriter` swallows a flush error by design. A file cut at an entry
  boundary is byte-for-byte the prefix of a longer one: the probe cut a three-entry spill to the
  end of its second and got **2 entries and 0 errors** where the writer had reported 3.
- **Why it matters:** pass three writes the VCF from the spill, so a spill silently short by its
  last entries produces a VCF silently short by its last records — on the one artifact whose
  byte-identity is the filter's standing oracle.
- **Fix:** correct the doc now, and carry the count from `entries_written` to the reader when B2
  owns the file's lifecycle.

- **M10: [spill.rs:95,98](../../../../src/ng/run/paralog_filter/spill.rs#L95-L98) — the two flags make a state the spec excludes representable, and nothing rejects it**
- **Categories:** idiomatic, refactor_safety — convergent
- **Confidence:** High on the shape; the consequence is a judgement about C3, which is unwritten.
- **Problem:** spec §3.2 puts "a repeat tract" inside "**every other record**", so
  tract-and-biallelic-SNP cannot occur — yet both flags are independent `bool`s with different
  consumers (the ordering check and the scorer), so a producer that set both would be caught at
  neither end. Separately, `is_biallelic_snp` sits on the entry and `alt_reads` on every sample
  with the coupling stated only in prose, and spec §6 trap 1 says C3 is where that gets confused:
  a scorer that reads `alt_reads` unconditionally compiles, runs, and returns a number.
- **Fix:** two agents proposed a three-variant sum type. **That changes spec §3.7's type sketch
  and is open question 1.** The half that needs no ruling is to refuse the impossible pair as
  corruption, alongside `NotABoolean`.

- **M11: [tests.rs](../../../../src/ng/run/paralog_filter/spill/tests.rs) — no property test over the codec, though `proptest` is already a dev-dependency**
- **Categories:** reliability
- **Confidence:** High.
- **Problem:** the reliability rules name serializers and round-tripping laws as requiring a
  property or fuzz test. All 17 tests are hand-written fixtures over four entry shapes. The
  classes the fixtures miss are exactly the ones a generator finds for free — a line at 127, 128
  and 129 bytes, a sample count at the same boundary, `-0.0`, both infinities. Two of those are
  M4, found by hand only because someone went looking.
- **Fix:** one property over generated entries. Worth adding only **after** B1, since a property
  test resting on a comparator that cannot fail proves nothing.

#### Minor

- **Mi1: [spill.rs:37-38](../../../../src/ng/run/paralog_filter/spill.rs#L37-L38)** — "Every field is
  either fixed-width or carries its own length" is false of the five varints, which are
  *self-delimiting*. This sentence is the module's whole justification for having no outer frame.
  (naming)
- **Mi2: [spill.rs:209-212](../../../../src/ng/run/paralog_filter/spill.rs#L209-L212)** — `finish`'s
  doc claims a guarantee the code does not provide; the wording holds in `vcf/writer.rs`, where
  `finish` performs the rename, and does not transfer. (naming; the behaviour is M9)
- **Mi3: [mod.rs:7](../../../../src/ng/run/paralog_filter/mod.rs#L7)** — "spills once" does not
  distinguish production from ng, which also spills once; per spec §2 the difference is spill
  *reads*, one against two. The sentence written to justify the module states something true of
  both sides. It also uses *spill* as a verb two lines before the file defines *the spill* as a
  noun. (naming)
- **Mi4: [spill.rs:66](../../../../src/ng/run/paralog_filter/spill.rs#L66)** — the doc cites
  `doc/devel/ng/spec/window_coverage.md` as the authority for where `WindowCoverage` belongs, and
  that file is not in the repository at this commit. (module_structure)
- **Mi5: [spill.rs:61-77](../../../../src/ng/run/paralog_filter/spill.rs#L61-L77) and
  [mod.rs:25](../../../../src/ng/run/paralog_filter/mod.rs#L25)** — the stand-in `WindowCoverage`
  is declared in the codec file and re-exported, so a type scheduled for deletion has two public
  paths, and the crate already has a `WindowCoverage` at `sample_summary/coverage.rs:274`.
  (module_structure)
- **Mi6: [spill.rs:337-341](../../../../src/ng/run/paralog_filter/spill.rs#L337-L341)** —
  `SAMPLE_CAPACITY_HINT` names a ceiling as a hint, justifies itself by pointing at a document
  rather than naming the number it beats (spec §4's 3,000), does not state its cost
  (`size_of::<SpilledSample>()` = 16 bytes, so 128 KiB at the cap), and is used nineteen lines
  before it is defined. (naming, defaults)
- **Mi7: [spill.rs:172-188](../../../../src/ng/run/paralog_filter/spill.rs#L172-L188)** —
  `SpillReader` states its buffering requirement at length and `SpillWriter` states none, yet
  `append` issues one `write_all` per entry: a bare `File` sink costs a write syscall per called
  record, and B2/C2 pick that sink against a constructor that never mentions it. (defaults)
- **Mi8: [spill.rs:360](../../../../src/ng/run/paralog_filter/spill.rs#L360)** — a ten-byte varint
  above `u64::MAX` is accepted with six bits dropped, where every sibling corruption is refused;
  the root cause is `psp/varint.rs:96`'s `data << 63`, which is out of scope, so the check belongs
  here. (reliability, errors, extras — convergent)
- **Mi9: [spill.rs:311,321](../../../../src/ng/run/paralog_filter/spill.rs#L311)** — two error
  labels are prose (`"line length"`, `"sample count"`) where the other nine are the field's own
  spelling; a third spelling, `"line"`, names the payload. (reliability, refactor_safety,
  idiomatic, smells — convergent)
- **Mi10: [spill.rs:312,322](../../../../src/ng/run/paralog_filter/spill.rs#L312)** — the reader
  allocates a fresh `Vec` for the line and one for the samples on every entry, while the writer
  keeps a reused scratch buffer and its field doc explains why. At 3,000 samples that is 48 KB
  per record, over two passes. (smells)
- **Mi11: [tests.rs:315-377](../../../../src/ng/run/paralog_filter/spill/tests.rs#L315-L377)** —
  `OutOfRange` on `ref_reads` and `alt_reads` is untested, so a `read_varint`-instead-of-`read_u32`
  slip in either would wrap silently. (extras)
- **Mi12: [tests.rs:225](../../../../src/ng/run/paralog_filter/spill/tests.rs#L225)** — the boundary
  test covers the top of every field but neither the bottom nor the floats' signed edges. `-0.0`
  is the one a bit-pattern codec can lose in a way `==` hides, which is the same argument the
  module already makes for `NaN`. (reliability)
- **Mi13: [tests.rs:174-196](../../../../src/ng/run/paralog_filter/spill/tests.rs#L174-L196)** — the
  "already carries a FILTER" and "INFO is a dot" round trips cannot take a path the base round
  trip does not: the codec writes the line with `extend_from_slice` and reads it back by length,
  so nothing branches on its contents. The plan names both cases, so keep them and say what they
  prove; the column-level versions are B3's. (smells)
- **Mi14: [spill.rs:84-86,96-97](../../../../src/ng/run/paralog_filter/spill.rs#L84-L97)** —
  `is_biallelic_snp`'s doc says the decision happens "here" and no code in this module makes it;
  and a record-level field's doc says "this sample's" where no sample is in scope. (naming)
- **Mi15: [spill.rs:1-9](../../../../src/ng/run/paralog_filter/spill.rs#L1-L9)** — *the verdict* and
  *the cut* carry the opening argument before either is defined, and "the two columns the verdict
  changes" never names them, where the spec does (`FILTER` and `INFO`). (naming)
- **Mi16: [spill.rs:11,242](../../../../src/ng/run/paralog_filter/spill.rs#L11)** — two clauses tell
  the reader a fact is important instead of stating it ("that is the design's load-bearing
  choice", "is the whole reason the distinction is drawn here"), which is a standing correction in
  `clear-technical-writing` Rule 6. (naming)
- **Mi17: [mod.rs:19,21](../../../../src/ng/run/paralog_filter/mod.rs#L19)** — "Milestone B" and
  "Milestone C" are plan-internal labels in the module's status line. (naming)
- **Mi18: [tests.rs:75](../../../../src/ng/run/paralog_filter/spill/tests.rs#L75)** — `same_bits` is
  a `-> bool` that does not read as a predicate, and its local `heads_match` covers six fields
  where `SpillEntry`'s own doc fixes *head fields* as the three the ordering check reads. (naming)
- **Mi19: [spill.rs:313](../../../../src/ng/run/paralog_filter/spill.rs#L313)** — `read` names a byte
  count fifteen lines below `ref_reads`/`alt_reads`, where *read* means a sequencing read. (naming)
- **Mi20: [spill.rs:244](../../../../src/ng/run/paralog_filter/spill.rs#L244)** — `next_entry` is a
  second public name for `Iterator::next` with the same body and no caller outside the module, and
  has no `# Errors` section where its siblings do. (idiomatic)
- **Mi21: [spill.rs:172](../../../../src/ng/run/paralog_filter/spill.rs#L172)** — `SpillWriter`
  carries a `finish` obligation nothing marks; `entries_written` has `#[must_use]` and the type
  does not. (idiomatic)

#### Nits

`gc_fraction`'s doc says "between 0 and 1" where the type doc says `NaN` is legal;
"leaves nobody anywhere to look" is ungrammatical; "the layout this module's header lays out"
repeats itself; the fixture helpers `present` and `without_a_window` are a bare adjective and a
modifier with no noun where every other fixture is a noun phrase; `encoded` is a past participle;
`TINY_TRACT_FLAG_AT` parses as "tiny tract flag"; `the_bytes_are_the_layout_the_spec_lays_out`
names no behaviour; the `NotABoolean` variant and the test that drives it use two words for one
concept (*flag byte* / *boolean*); `mod.rs`'s title names a shape rather than a domain thing;
`SAMPLE_CAPACITY_HINT` is declared after its only use; `as u64` on the encode side where the
decode side uses `try_from`; `assert_eq!(read.len(), 3)` then three index expressions where a
slice pattern states the shape once; and no `TINY_BIALLELIC_FLAG_AT` beside `TINY_TRACT_FLAG_AT`,
which is the mechanical reason that flag is the one no test corrupts.

### 7. Out of scope observations

- **[varint.rs:93-112](../../../../src/psp/varint.rs#L93-L112)** — `decode_u64_leb128_cold` folds
  the tenth byte in as `data << 63`, so a ten-byte varint above `u64::MAX` returns `Ok` with six
  bits discarded rather than `VarintError::Overflow`, and it also accepts a non-canonical
  zero-padded encoding of a small value. The `.psp` reader is its other caller. Production's
  module; the guard for this plan's use goes in `read_varint` (Mi8).
- **Spec §3.4's heading** still reads "The spill — a **framed** binary file of finished lines"
  while the implementation has no per-entry frame. The report records the deviation and argues
  it; the spec's own heading says the opposite, which is a spec edit or a ruling, not a code fix.

### 8. Missing tests to add now

Grouped by what they defend. Bodies are in the per-category files.

**`same_bits`** — `same_bits_separates_entries_that_differ_in_any_one_field`: two entries
differing in exactly one field, one pair per field, including an absent window against a zeroed
one. Catches a comparator that has stopped discriminating (B1).

**`decode_entry`** — `a_file_cut_anywhere_inside_a_record_names_the_field_the_bytes_ran_out_in`:
every truncation offset of a known encoding with its expected label (M5).
`a_biallelic_flag_byte_that_is_neither_zero_nor_one_is_refused` (M5).
`a_line_longer_than_the_ceiling_is_refused_before_it_is_read` (M1).
`a_ten_byte_varint_past_u64_is_refused_rather_than_truncated` (Mi8).
`a_read_count_too_large_for_its_field_is_refused` (Mi11).

**The round trip** — `a_line_longer_than_one_varint_byte_round_trips` and
`a_cohort_of_a_thousand_samples_round_trips` (M4).
`any_stream_of_entries_comes_back_bit_for_bit`, a `proptest` property (M11).
The boundary fixture extended with `-0.0` and both infinities (Mi12).

**`SpillReader`** — `the_reader_stops_after_a_decode_error`: a corrupt prefix followed by a
well-formed entry (M3).

**`SpillWriter`** — `append_reports_a_sink_that_refused_the_bytes_and_does_not_count_the_entry`
and `finish_reports_a_flush_that_failed`, over a sink that refuses everything (M7).

### 9. What's good

- **The pinned byte list** (`the_bytes_are_the_layout_the_spec_lays_out`) is the only test that
  catches the two floats swapped on *both* sides; the nine round trips are blind to it by
  construction.
- **The encoder's exhaustive destructure works as advertised** — a third field on
  `WindowCoverage`, `#[non_exhaustive]` on top of it, and a rename each produce compile errors at
  the encoder, the decoder and the fixtures. The stand-in cannot lose a field at the rebase.
- **The `SAMPLE_CAPACITY_HINT` clamp is the right shape**: it changes allocation only, never a
  decoded value, so a cohort above it decodes identically and pays one `Vec` growth. Removing it
  aborts the process on the four-billion fixture with
  `memory allocation of 64424509440 bytes failed`.
- **The `BufRead` bound is not a preference but the only stream-safe way** to use a slice-taking
  varint primitive across a buffer boundary, and the writer's reused scratch buffer is forced by
  `encode_u64_leb128` taking a `&mut Vec<u8>`.
- **The module's decisions are argued where they are made**, and three of them were checked
  against the spec and found faithful — the layout field for field, the float widths, and the
  absence of a per-entry frame.

### 10. Commands to re-verify

Run in the container from the branch's own `scripts/dev.sh`:

- `cargo test --all-features --lib "ng::run::paralog_filter::"`
- `cargo test --all-features --lib --bins --tests`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings`
- `cargo fmt --check` — and `git checkout --` the nine main-owned files afterwards

New, introduced by this review: the nine tests of section 8.
