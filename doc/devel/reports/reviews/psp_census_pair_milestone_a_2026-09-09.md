# Code Review: the census lives inside the psp — Milestone A
**Date:** 2026-09-09
**Reviewer:** rust-code-review skill (orchestrator), two category agents
**Scope:** the working-tree diff of each step of Milestone A, one review per step
**Status:** Approve-with-changes — every finding below was applied in the step that raised it

---

## How this review was run, and where it departs from the skill

**Two read-only agents in the main checkout, not a per-category fan-out in worktrees.** The skill
asks for one agent per category, each in its own `isolation: "worktree"` tree, so they can mutate
freely. Measured on this machine: the project's build tree is **17 GB** and the disk has **129 GB
free**, so a seven-agent fan-out with a cold build each would not fit, and four steps of it would
not fit twice over. The diff under review is two files.

So the categories were grouped into two agents — *reliability + errors*, and
*unsafe_concurrency + smells + refactor_safety* — both forbidden to edit any file or run `cargo`,
and both asked to name any mutation they wanted rather than run it. **The orchestrator ran the one
mutation they asked for, serially, and it changed the outcome** (see A1's finding R3). That is the
trade this arrangement makes: fewer mutations, each run by one hand, against a fan-out that could
have run more.

---

## A1 — the walk hands the census to `finish`

**Reviewed against:** the working tree over `9fdc427f`, files `src/ng/run/psp_writer_line.rs` and
`src/ng/run/gatherer.rs`.

### Verification given to the reviewers

`cargo check --lib --bins --tests` clean but for the baseline's own warning;
`cargo test --lib psp_writer_line` 5 passed; `cargo test --lib gatherer::census_tests` 6 passed;
`cargo fmt --check` flagging the baseline's four files and neither of these two. The tree's
pre-existing red is recorded in the
[implementation report](../implementations/ng_psp_census_pair_milestone_a_2026-09-09.md).

### Findings

**M1 — `src/ng/run/gatherer.rs:461` — `write_psp`'s own doc still said the trailer is sealed empty.**
The first sentence of the function this change is about read "the trailer is sealed **empty**, a
legal trailer whose contents stay opaque until something needs them", citing the container spec.
That is the belief the design exists to end, in the doc a reader meets first.
*Fixed:* the summary now says the seal carries the sample's census, and cites `psp_census_pair.md`
§3.1.

**M2 — `src/ng/run/psp_writer_line.rs:182` — `Vec::new()` and a forgotten census were the same value
to the compiler.** `finish(trailer: Vec<u8>)` gives no way to tell *this walk built no census* from
*the caller forgot one*, and `line.finish(Vec::new())` was already the idiom at three of four call
sites. The product is a psp that opens, indexes and reports a correct record count with no census
in it — the state spec §3.3 exists to remove, found at the parameters fit.
*Fixed:* a `PspTrailer` enum with `Census(Vec<u8>)` and `Nothing`, local to the module, so the psp
format layer is untouched and closing with nothing is a sentence the author writes.

**M3 — `src/ng/run/gatherer.rs:540` — `CensusNotWritten` was the wrong error for an in-memory
encode.** It renders as *the census at …/zeta.psp could not be written*: the operation was an
encode into memory, the path names a different file, and at that moment the file at that path is a
half-written psp.
*Fixed:* a new `RunError::CensusNotEncoded { sample, source }`, naming the sample and no path. The
misfiled doc block above `CensusNotPlanned` — it described `CensusNotWritten` and ran into its
neighbour's own first line — moved onto the variant it is about.

**M4 — `src/ng/run/gatherer.rs:1650` — the trailer test could not detect one of the three failures
its doc claimed.** Both sides of its equality are encodings of the *same* `SampleCensusEvidence`,
so a census closed before the walk fed it is equally empty on both sides and compares equal. The
reviewer named the mutation: move the `CensusWriter::finish` above the walk loop.
*Ran it.* `gatherer::census_tests` stayed at 6 passed; `census_from_psp`'s two parity tests went
red, which is what the reviewer predicted and what proves the mutation bites.
*Fixed:* the overclaim removed from the doc, and an assertion added that some kept position has a
read at it. **The first version of that assertion did not work either** — it tested *not never
walked*, and building the writer marks every kept position of the analysed ground as walked
(`gatherer.rs:294`), so an unfed census is walked everywhere at depth zero. Rewritten as *depth
above zero*, the mutation fails the test. Reverted, re-run, restored tree confirmed by reading the
`+`/`-` lines of `git diff HEAD -- src/ng/run/gatherer.rs`.

**M5 — no test wrote a trailer over 31 bytes, and this step makes production trailers tens of
megabytes.** Every `PspWriter::finish` under test passed an empty slice or a short literal. A
truncation or offset defect above a buffer boundary would leave a psp that opens, indexes and reads
back every record and hands the fit a census that stops early.
*Fixed:* `a_trailer_of_megabytes_round_trips_through_the_line`, 4,000,000 varying bytes.

**M6 — a failure writing the census beside the psp now discards a psp that is whole.** By the time
`write_census_beside` runs, the census is already in the trailer — but `generate-psps` deletes the
part-written psp on any error out of `write_psp`, which was right while a psp without a census was
useless. It cannot be fixed by reordering: the sidecar's identity needs `stats.header_digest`,
which exists only once `finish` has returned.
*Recorded, not fixed.* Teaching `generate-psps` the difference is work in a file this step does not
touch, for a failure mode that exists only until A2 deletes the copy. The ⚠ paragraph on
`write_census_beside` says so.

**Mi1 — `write_psp`'s `# Errors` omitted the error it can now raise two ways**, and the two leave
opposite things behind: an encode failure leaves no readable psp, a sidecar failure leaves the psp
whole and carrying its census. *Fixed.*

**Mi2 — `write_loci`'s seal is a variable assignment and the loop carries on reading.** Unreachable
today — `finish` is the only sender, it consumes the line, the sender is never cloned — but the
shape allows a second seal to replace the first payload silently, and a record arriving after a
seal to be written and sealed over. *Fixed* with two `debug_assert!`s naming both invariants.

**Mi3 — the discarded `send` result now carries the only copy of the census**, with nothing at the
site saying why the discard is safe. *Fixed:* the one-line reason — a failed send means the thread
panicked, and the `join` below re-raises it.

**Mi4 — `census_as_a_trailer` documented an error no test can cover.** A `Vec<u8>` as an
`io::Write` never fails. *Fixed:* the doc now says so, so the next reader does not go looking for
the missing test.

**Mi5 — the trailer test's oracle is the file A2 deletes.** *Recorded on the test*, naming A4's
rebuild-from-the-psp as its replacement.

**Nits, all fixed:** the "encoding into memory can only fail the way a `Vec` fails" shorthand, which
states its conclusion in terms the reader has to unpack; "Until this line" read as *up to this line
of code* where it meant *before this change*; the module note's "crosses the seam once", which a
reader would take to mean the census is encoded once, when in this step it is encoded twice; a test
naming its file `alpha.psp` while walking zeta; and a `let _ =` binding a value `expect` had already
unwrapped.

### Carried forward, not fixed

- **`write_census_beside`'s two independent `Option`s**, whose fourth combination cannot happen and
  is caught by a runtime `assert!`. The function is deleted in A2.
- **`write_loci`'s `Ok(Ok(None)) => unreachable!`** — a panic standing in for an invariant.
  Pre-existing and still sound.
- **`RunError::CensusNotWritten`'s `Box<dyn Error>` source** — pre-existing; the new
  `CensusNotEncoded` beside it carries a typed source instead.
- **Memory, recorded as a number rather than a fix.** `write_census` materialises every section
  before writing a byte, so encoding the trailer holds the evidence, the sections and the
  concatenated bytes at once, and the bytes stay live across the seal. At spec §3.1's whole-genome
  size — 1.03 bytes a position and about 35 a tract — that is a transient of tens of megabytes, and
  this step's sidecar encodes it a second time. A2 removes the second encode; the trailer's copy
  stays until `write_census` learns to stream into the psp writer.
- **`reader.rs:236` allocates the trailer in one go from a length the file supplies.** Harmless
  while every trailer is empty; A1 makes the field large by design and A3 reads censuses out of psps
  this process did not write. It belongs to A3's review.

---

## A2 — `generate-psps` writes one file

**Reviewed against:** the working tree over `e1ae9d98`, 13 files. Same arrangement as A1: two
read-only agents, grouped categories, mutations named rather than run. **Three were named and all
three were run serially** — see the implementation report's table; each was reverted and the
module's tests re-run on the restored tree.

### Findings

**M1 — `generate_psps.rs` — the read-back after the rename was reported as a stopped walk.**
`census_bytes` was read by reopening the finished psp, and a failure of that open became
`GeneratePspsCliError::Walk`, whose own documentation says the sample has no psp and that the psp
it was replacing is untouched. On that path the walk had finished and the rename had already
happened, so both were false, and the message told the reader to re-walk a sample that had
succeeded.
*Fixed differently from either suggestion, and it dissolves the finding:* `WriteStats` carries
`trailer_bytes`, so the writer hands back the number it already had. The reopen, its error path,
and the report's two numbers coming from two sources all go together.

**M2 — the report's census clause was asserted nowhere.** Nothing in the tree looked for the
strings the per-sample line and the run's totals print; `the_per_sample_line_prints_the_numbers_it_names`
set `census_bytes: 812` and never looked for 812, and the report test compared the struct field
against the same `footer().trailer_bytes` call the production code had used.
*Fixed:* both strings are asserted, the trailer is read through `trailer()` and decoded rather than
re-derived from the footer, and the totals line is asserted with a fixture check that the two
totals differ. Two mutations that were green are now red.

**M3 — the record count in a census file's identity was covered nowhere** once the command-level
test that counted a psp's records went with the sidecar. The cohort opener compares only the header
digest.
*Fixed:* `the_census_it_writes_names_the_psp_it_read`. A `records: 0` mutation was green and is now
red.

**M4 — `each_census_it_writes_equals_the_one_the_walk_wrote` was weakened without saying so.** It
decodes the file `generate-census` wrote and re-encodes it without its pileup identity, so the
file's own layout survives only as far as `decode_census` preserves it, the identity is dropped
unchecked, and the test newly leans on a byte-level round-trip law that `census_file.rs`'s tests do
not pin — they compare decoded values.
*Fixed:* the doc says what is compared; the round-trip law gets
`write_census_after_decode_census_returns_the_bytes_it_was_given`; the identity gets M3's test.

**M5 — two shipped scripts consume the file this step deletes.** `scripts/ng_census_route_cost.sh`
and `scripts/ng_fit_stage_end_to_end.sh` both glob `*.census` and both now fail on every run, the
second before it reaches step 3.
*Recorded, not fixed.* The plan puts them in step E1, whose other half needs commands that do not
exist yet. **The consequence is that the end-to-end harness cannot verify anything on real reads
between here and Milestone E**, and that is stated in the implementation report rather than left
for whoever runs it next to discover.

**Minor, all fixed:** `a_stopped_walk_leaves_neither_file_at_the_samples_own_path` named a pair that
no longer exists and never asserted the guarantee its doc stated — renamed, and it now walks the
cohort once first; `force_replaces_a_census_that_is_already_there` compares `generate-census` with
itself while its doc explained it as two producers agreeing — reworded, with a pointer to where the
cross-producer guarantee now lives; five doc claims in `gatherer.rs` this change falsified, the
module doc of `generate_psps.rs`, and `examples/ng_census_route_cost.rs`'s header, which said the
two routes produce the same file beside the psp; a discarded `remove_file` result with no reason
given; and the totals line's "bytes of it census", whose *it* had no referent and which phrased the
same relation differently from the per-sample line.

**`CENSUS_FILE_EXTENSION` and `census_path_for` in `generate_psps.rs` are left where they are.**
Both are now test-only and their docs say "this command writes"; plan step E2 deletes them, and
moving them a step early would touch four modules for no behaviour.

**The fixture helper was reshaped rather than accepted.** `censuses_written_beside_the_psps` first
took a reference, a catalog and a directory, and hardcoded the five repeat criteria as defaults —
which every caller happens to use, so it would have agreed with all of them until one walked with
anything else. It takes the `GeneratePspsArgs` the psps were walked under. `build_every_census`,
widened to `pub(crate)` for it, is private again: `run_generate_census` was already public and does
the same.

---

## A3 — the fit reads a census out of the middle of a file

**Reviewed against:** the working tree over `1906752c`, three files. One read-only agent over five
grouped categories, given the mutation I had already run so it would not spend a round re-finding
it.

### Findings

**M1 — the doc's account of what the length buys was wrong in direction and in size.** It said an
uncapped head read "would pull a megabyte of records past each census's end". In a psp the records
come *before* the trailer; `read_footer` refuses any file whose trailer does not end exactly where
the footer begins, and the footer is 48 bytes — so an uncapped read runs at most 48 bytes past and
then hits the end of the file. And at the size the same paragraph quoted, the cap is a minimum
against a 1 MiB buffer and does not bite at all.
**This was the second wrong account of that argument in one step** — the first, corrected before
the review, said the cap prevented a malformed decode, which it does not either.
*Fixed by giving the length a job it actually does* (M2), and by saying plainly that the head-read
cap is a size choice with no test.

**M2 — `len: Option<u64>` let the caller drop the argument the function exists for, silently.**
Every caller that has the offset has the length beside it in the same footer, and the module
already owns the type for that pair.
*Fixed:* `open_census_within(path, census: ByteExtent)`. And the length now bounds the directory at
open: **every section must end inside the census**, so a directory that outgrew it is refused
rather than turning into a seek into the psp's records and a `resize` to a length the file
supplied — which for a large enough value aborts the process instead of returning an error.
`a_section_that_ends_outside_the_census_is_refused_at_the_door` is what fails without it.

**M3 — nothing tested that relabelling a census carries its offset**, and relabelling is on the
path every sample of every cohort takes. Both of the step's tests read an un-relabelled census, and
every other cohort test opens at offset zero.
*Ran the mutation the reviewer named* — writing `at: 0` in that arm — *and it was green.* Now
`a_census_in_a_cohort_keeps_where_it_is_in_its_file` fails on it: three samples, each at a
different offset, read through the cohort after it has renumbered them.

**Minor — what is still not bounded, stated rather than fixed.** A section's extent is now checked
against the census's length; a section's *content* is not, and what catches a misdirected read is
`decode_section`'s exact-consumption rule rather than the directory's overlap check, which compares
sections against each other and cannot see the census's end. On files this project wrote the live
risk is corruption, not attack, and the outcome is a refusal. Recorded in the implementation
report.

**Minor, fixed:** the module header said a census lives "beside that sample's pileup and never
inside it" and that offsets are "from the start of the file"; `Sections` said "only the resident
state exists today"; and the byte-count assertion claimed a section's read excluded the padding
around the census, when the counter does not cover the head read at all.

**Confirmed by the reviewer and not changed:** every match on `Sections::Backed` is correctly
updated, and `Sections::backed` is the only other construction site.
