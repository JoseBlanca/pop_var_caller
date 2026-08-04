# Code Review: ng_read_filtering_stages_a1
**Date:** 2026-08-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A1 of `impl_plan/read_filtering_stages.md` — `RawRecord` → `RawAlignedRead`,
`NoodlesRawRecord` → `NoodlesRawAlignedRead`, both moved into `read/aligned_read.rs`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted working-tree diff for plan step **A1**, exported as
  `tmp/review_2026-08-03_ng-read-filtering-stages-a1/a1.patch` and re-applied by each agent onto
  a detached `8cf6f03`.
- **Reviewed against:** branch `ng-generic-perf`, base commit `8cf6f03`.
- **In-scope files:** `src/ng/read/aligned_read.rs`, `src/ng/read/filtering.rs`,
  `src/ng/read/mod.rs`, `src/ng/read/input/mod.rs`, `src/ng/read/input/region_records.rs`,
  `src/ng/read/input/record_reader/{mod,bam,cram,in_memory}.rs`,
  `doc/devel/reports/implementations/ng_read_filtering_stages_a1_2026-08-03.md`, and the new
  `PROJECT_STATUS.md` block.
- **Deliberately out of scope:** steps A2/A3 (the reader and region renames) and Milestones B–D
  (the contig check, the loop move, the deletions) — later commits in the same plan;
  `read/filtering/` as a folder (spec §10 defers it).
- **Categories dispatched:** `naming` (the step is a rename), `refactor_safety` (the step is a
  move and must change nothing), `module_structure` (the step re-homes a type),
  `reliability` (what the moved code does and does not pin), `extras` (diff-matches-intent, and
  the impl report's accuracy).

## 2. Verdict

**Approve-with-changes.** The rename and the move are clean and provably behaviour-free. What
the review found is not drift but **absence**: three of the invariants the moved code documents
at length are pinned by nothing, and the test that looks like it guards them passes for the
wrong reason.

## 3. Execution status

Run by the orchestrator on the main checkout (host, debug) and independently reproduced by four
of the five agents in their own worktrees:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean, 8.56s |
| `cargo test --lib` | 0 | 2,837 passed / 0 failed / 5 ignored — **unchanged from base** |
| `cargo test --lib ng::` | 0 | 1,538 passed / 0 failed / 2 ignored — **unchanged from base** |
| `cargo test --examples` | 0 | 52 passed / 0 failed |
| `cargo doc --no-deps` | 1 | 12 unresolved links, **all pre-existing, none in a touched file** |

**Not run, with reasons already on record:** `cargo test --release` (red on a clean tree — four
tests assert on `debug_assert!` messages release compiles out); `cargo test --all-targets`
(aborts on a pre-existing panic in `benches/psp_writer_perf.rs:386`); `cargo audit` (not part of
this project's gate).

**The four acceptance dumps are byte-identical to the `8cf6f03` baseline** by `cmp` — 251,792 /
4,406 / 1,718,914 / 11,945 lines — and `ng_generic_walk_probe` prints the anchor exactly:
`loci=236081 observations=251786 reads_admitted=54709`.

Findings labeled "Needs verification": **0.** Every finding below was produced by running a
command, not by reading.

## 4. Open questions and assumptions

1. **Should `NoodlesRawAlignedRead` be narrowed from `pub` to `pub(crate)`?** (affects Mi2)
   The `module_structure` agent verified the narrowing plus dropping the dead re-export leaves
   `clippy -D warnings` clean. It is a **crate public-API change**, not a rename, so it is not
   this step's to take silently — raised for the owner at Checkpoint A.
2. **Were the two Blocker fixes in scope for a renames-only milestone?** Both gaps are
   pre-existing: the moved impl is byte-identical to `8cf6f03`. They are filed here because A1 is
   the commit that re-homes this code beside the tests meant to guard it, and because a
   renames-only milestone is exactly when you find out what was never tested. The alternative
   offered by the `refactor_safety` agent — park them as a rider on C2 — is recorded, not taken.
   (affects B1, B2)

## 5. Top 3 priorities

1. **B1** — `decode`'s refusal of an unstamped buffer can be replaced by the exact silent
   fallback it exists to prevent, and all 2,837 tests pass.
2. **B2** — the hand-written `Default` can be given the trap value its own eight-line doc
   forbids, and all 2,837 tests pass; the trait accessor is equally unpinned.
3. **M1** — the one test that looks like it covers the above passes for the wrong reason: it
   never reaches the path its name, its comment and the impl report all claim.

## 6. Findings

### Blocker

**B1: src/ng/read/aligned_read.rs:255-266 — `decode`'s read-group guard has no test that fails when the refusal is deleted**
**Categories:** reliability · **Confidence:** High

`decode` refuses a buffer whose `read_group` is `None`, because `ReadGroupId(0)` is a real id and
a missed stamp would otherwise attribute every read to a real library. The guard was replaced
with the exact silent fallback it exists to prevent —
`let read_group = self.read_group.unwrap_or(ReadGroupId(0));` — and `cargo test --lib` returned
`2837 passed; 0 failed`. Nothing in the crate notices.

The test that looks like it covers this passes on *both* branches: with the guard it fails on
the missing stamp; without it, the default `RecordBuf` fails on its missing reference sequence
id — and the test asserts only `err.kind() == InvalidData`, which both produce.

**Failure scenario:** a reader arm added later forgets to stamp the buffer. Every one of its
reads is attributed to `ReadGroupId(0)`, a genuine library. Per-read-group drop tallies and
everything fitted per read group are wrong, with no error and no dump difference.

**Fix:** a test whose fixture is *decodable*, so the missing stamp is the only thing that can
fail it. Verified by the agent: passes unmutated, fails under the mutation.

**B2: src/ng/read/aligned_read.rs:223-238, 251-253 — the hand-written `Default` and the `read_group()` accessor are both unpinned**
**Categories:** reliability, refactor_safety (convergent) · **Confidence:** High

The `Default` impl carries the longest justification in the moved code: "`ReadGroupId(0)` would
be the obvious filler and is the trap: zero is not a sentinel, it is the first identifier the
run's table mints … `None` makes that a loud failure at the one place that reads the field."
Both agents independently wrote that trap in — `read_group: Some(ReadGroupId(0))` — and both got
`2837 passed; 0 failed`. Separately the accessor, mutated to
`Some(self.read_group.unwrap_or(ReadGroupId(0)))` so it can never report an unstamped buffer,
also survives all 2,837.

The existing `in_memory.rs` and `region_records.rs` tests assert on the **field**
(`buf.read_group.is_none()`), never through `RawAlignedRead::read_group()`, and none constructs
a buffer via `Default` and inspects it.

**Failure scenario:** with `Default` stamping `Some(ReadGroupId(0))` the guard at `decode` can
never fire for a fresh buffer, so the whole `Option` becomes decorative — B1's failure, reached
from the other side.

**Fix:** one test that reads a defaulted buffer *through the trait accessor*, killing both
mutations at once. Verified by the agent against each.

### Major

**M1: src/ng/read/aligned_read.rs:562-569 — `..._decode_errors_on_a_record_with_no_position` tests neither thing its name and doc claim**
**Categories:** reliability, refactor_safety (convergent) · **Confidence:** High

The test drives `NoodlesRawAlignedRead::default()`, whose `read_group` is `None`, so `decode`
returns at the read-group guard and `decode_record` is never called. Its comment says the
opposite: "A default record has no reference_sequence_id / alignment_start, so the reused
decoder fails". Both agents proved the real path — one by probe
(`the error this test actually sees is: a record reached decode with no read group`), one by
mutating only the guard's error kind and watching this test fail on `NotFound` vs `InvalidData`.

So the adapter's *actual* corrupt-record path — a stamped buffer holding a record with no
alignment start — is exercised by nothing, while the test's name says it is. This is also what
let B1 and B2 survive.

**Failure scenario:** a future edit deletes the wrong one of the two tests, believing the
corrupt-record path is covered here rather than at
`a_record_with_no_position_fails_naming_the_read`.

**Fix:** repurpose it to the subject its name promises — stamp the buffer so the decoder's own
error is the one that surfaces — and assert on the message, not just the kind.

### Minor

**Mi1: src/ng/read/input/mod.rs:26-29 — the one prose edit outside the moved files now misattributes the type**
**Categories:** naming, module_structure, extras (three-way convergent) · **Confidence:** High
The doc reads "it lives under `read/` beside `filtering.rs` and reuses **that module's**
`RecordSource`/`RawAlignedRead` seam". After A1 only `RecordSource` is `filtering.rs`'s. The
diff renamed the link target but not the sentence around it, and `cargo doc` stays quiet because
`super::RawAlignedRead` still resolves through the re-export. It is the one place in the diff
that still says the types did not move.

**Mi2: src/ng/read/aligned_read.rs:207, src/ng/read/mod.rs:35 — the move carried `pub` across unexamined**
**Categories:** module_structure · **Confidence:** High
`NoodlesRawAlignedRead` is `pub` with no consumer outside `src/ng/read/`, and the
`read::NoodlesRawAlignedRead` re-export has no consumer at all. The agent verified narrowing to
`pub(crate)` and dropping the re-export leaves clippy clean. It also verified the trait
`RawAlignedRead` **cannot** be narrowed — rustc: *"trait `RawAlignedRead` is more private than
the item `RecordSource::Record`"* — its `pub` is held open by `pub trait RecordSource`, a
Milestone C deletion. **Open question 1**: this is a public-API change, not a rename.

**Mi3: src/ng/read/aligned_read.rs:190, 202 — the moved type doc-links back into the module it just left, and mislabels whose trait it is**
**Categories:** module_structure · **Confidence:** High
An intra-doc link `[`RecordSource::read_next`](crate::ng::read::filtering::RecordSource::read_next)`
from the data-model module into the pipeline-stage module — the one back-reference the move
leaves, and a scheduled `cargo doc` breakage when Milestone C deletes `RecordSource`. The same
block also says "The **production** `RecordSource`"; `RecordSource` is **ng's** trait, which is
the exact fact the next paragraph exists to pin down.

**Mi4: src/ng/read/input/record_reader/in_memory.rs:206, bam.rs:374 — `..Default::default()` in struct literals whose subject is the field being spread past**
**Categories:** refactor_safety · **Confidence:** High
The checklist's explicit rule. `NoodlesRawAlignedRead` has two fields, so the spread fills
exactly one, saving nothing while removing the compiler's ability to flag these sites if the
struct gains a field. Both literals exist to plant a stale `read_group` and prove the reader
clears it — the spread is precisely the construct that would silently absorb a third field, and
Milestone C moves this buffer onto `AlignmentCursor`.

**Mi5: src/ng/read/aligned_read.rs:508-523 — `record_with`'s name and quality scores are inert**
**Categories:** reliability · **Confidence:** High
Both parameters are load-bearing (mutating either kills the consuming test), but
`.set_name(b"r1")` and `.set_quality_scores(vec![30u8; 4])` are asserted by nothing — changing
them leaves `ng::read::aligned_read` at `7 passed`. `decode`'s pass-through of `qname` and
`qual` *through the adapter* is unasserted.

**Mi6: the impl report is wrong in four checkable places**
**Categories:** refactor_safety, extras, naming (three-way convergent) · **Confidence:** High
(a) "**Five** call-site files" followed by a list of six. (b) `bam_record` "was the helper for
exactly those three tests … `clippy -D warnings` said so" — it had **one** call site, and clippy
reports dead code only after the last user leaves, so it never said that. (c) the trait is
"unchanged apart from its name and the new unmapped-read paragraph" — its **opening sentence was
rewritten too**, to arch §3.1's wording. (d) "its **six** now-unused noodles imports" — five
`use` lines carrying seven items. On a milestone whose value proposition is *the diff is exactly
the rename*, the report is what a reviewer reads instead of re-deriving the diff.

**Mi7: src/ng/read/aligned_read.rs:95 vs :62 — mapping quality is `MapQual` on the trait and bare `u8` on the struct, now 30 lines apart**
**Categories:** naming · **Confidence:** High
The same "SAM `0xFF` → 0" rule is now written twice at two result types. **Explicitly not an A1
edit** — `AlignedRead.mapq` is a public field and every reader of it would change. For B or C,
whichever first edits the struct.

**Mi8: src/ng/read/aligned_read.rs:42 — pre-existing back-reference from ng's read data model into production's `pileup` stage**
**Categories:** module_structure · **Confidence:** High
`use crate::pileup::walker::CigarOp;` — a stage-internal type consumed by ~19 files under
`src/ng/` alone. Not introduced by this diff and fixing it touches frozen production code.
For `module_layout.md`'s Open items, not for this commit.

### Nits

Collected, not enumerated: two consecutive `use` lines from the same module in
`region_records.rs:46-47` where the diff split one braced import; a now-redundant
`use noodles_sam::alignment::RecordBuf;` in `aligned_read.rs`'s test module (the diff's new
top-level import reaches it through `use super::*`, and rustc does not warn on
explicit-shadows-glob); a `FLAG_PAIRED` shadow the diff created (the new import vs two
function-local `const FLAG_PAIRED` in `readthrough_pair`/`readthrough_pair_reverse`, plus
`FLAG_MATE_REVERSE`/`FLAG_REVERSE` as second names for production's `*_STRAND` constants);
`err` where the merged test module uses `error`; `record_with`'s dangling name and its plural
"these tests" against one caller; an ungrammatical sentence in the new module doc ("What decides
*which* reads get converted is `read/filtering.rs`'s, and it stays there"); the now half-false
standing comment in `filtering.rs` ("This module knows what a *record* is"); two comment lines
the substitution pushed past their file's wrap (`read/mod.rs:48` at 106 columns,
`input/mod.rs:28` at 88); `FakeRecord` keeping the retired vocabulary while its doc took the new
one (moot if C3 deletes it); bare `raw` bindings where `raw_read` carries the noun.

## 7. Out of scope observations

- `doc/devel/ng/README.md:106` still describes `filtering.rs` as holding "the cascade, the
  `RecordSource`/`RawRecord` seam, the `ReadFilter` iterator" — stale in both the name and the
  home. Worth a sweep at Checkpoint A, when all three renames have landed.
- `doc/devel/ng/arch/read_filtering.md`, `spec/read_filtering.md`, `arch/alignment_file.md` and
  `arch/read_groups.md` still use `RawRecord`/`NoodlesRawRecord` and cite `filtering.rs` line
  numbers A1 invalidated. They record superseded designs and the new spec's §6 table is the
  migration record, so a pointer at the top of each may be all that is wanted — owner's call.
- `src/ng/read/input/region_records.rs:588` — `a_yielded_record_carries_its_read_group` asserts
  only `buf.read_group.is_some()`, never *which* id was resolved. The third weak link in the same
  read-group chain; worth a distinctive-id assertion when A2 rewrites this file.
- `src/ng/read/input/cursor.rs` and `filtering.rs` accumulate tombstone comments marking where
  deleted code stood — three in `filtering.rs` now. Milestone C deletes most of what surrounds
  them; worth one sweep then.

## 8. Missing tests to add now

Grouped by function under test. Each was written and run by the `reliability` agent in its own
worktree, and each was verified to **fail** under the mutation it names.

**`NoodlesRawAlignedRead::decode`**
- `noodles_raw_aligned_read_decode_refuses_a_buffer_with_no_read_group` — input class: a buffer
  valid in every respect *except* the stamp (the existing test's buffer is invalid in several at
  once, which is why it cannot discriminate). Catches: deletion or weakening of the read-group
  guard. Asserts the fixture is decodable first, or it proves nothing.
- `noodles_raw_aligned_read_decode_errors_on_a_stamped_record_with_no_position` — input class: a
  *stamped* buffer holding a corrupt record, the combination no test currently builds. Catches
  `decode` swallowing or mis-reporting the decoder's error, and restores the path M1's test name
  promises.

**`NoodlesRawAlignedRead::default` / `RawAlignedRead::read_group`**
- `noodles_raw_aligned_read_default_reports_no_read_group` — input class: a freshly
  `Default`-constructed buffer, read *through the trait accessor* (which is what the tally
  calls, and what the field-level assertions elsewhere miss). Catches both B2 mutations.

**Filed but not taken in A1** (a boundary test, not a Blocker fix; recorded so it is not lost):
`noodles_raw_aligned_read_reads_the_maximum_available_mapping_quality` — 254 is the largest
value that is a quality, and `MappingQuality::new(255)` is `None`, which is *why* `0xFF` never
arrives as a number. Both sides of the top of the range are currently untested.

## 9. What's good

- **The move is provably free.** The `refactor_safety` agent extracted the four moved blocks from
  the base file, applied `sed 's/RawRecord/RawAlignedRead/g'` and diffed: the only hunk is a
  three-line doc re-wrap. That is the standard a renames-only milestone should be held to.
- **The new module doc opens on the right idea** — one read in two states — and says plainly that
  no keep-or-drop rule lives there, which is the rule the rest of the plan depends on.
- **The unmapped-read paragraph gives the reader the reason not to "fix" the name** (SAM calls
  every line an alignment record), rather than just asserting the name is right.
- **Every call site reaches the moved type by its crate-absolute defining path**, not through the
  `read::` re-export and not by a `super::super::` walk.
- **The old names are gone from Rust entirely** — `grep -rn 'NoodlesRawRecord\|\bRawRecord\b' src
  tests benches examples` returns nothing, including string literals and test names.

## 10. Commands to re-verify

```
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --lib
cargo test --lib ng::
cargo test --examples
cargo test --lib ng::read::aligned_read      # the moved tests and their new siblings
```
Plus the four acceptance dumps and the walk probe, compared with `cmp` against the `8cf6f03`
baseline — the only check that can see a behaviour change this suite cannot.
