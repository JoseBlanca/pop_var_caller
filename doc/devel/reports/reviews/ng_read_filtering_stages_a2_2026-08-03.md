# Code Review: ng_read_filtering_stages_a2
**Date:** 2026-08-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A2 — `RecordReader` → `AlignedReadsReader`, its three arms, and the module
`record_reader/` → `aligned_reads_reader/`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted working-tree diff for plan step **A2**, exported as
  `tmp/review_2026-08-03_ng-read-filtering-stages-a2/a2.patch` and re-applied by each agent onto
  a detached `5438927` (A1).
- **In-scope files:** `src/ng/read/input/aligned_reads_reader/{mod,bam,cram,in_memory,container}.rs`,
  `src/ng/read/input/{mod,cursor,region_records,open_bam,sample_cursor,test_fixtures}.rs`,
  `src/ng/read/filtering.rs`, and the impl report.
- **Deliberately out of scope:** step A3 (the next commit) and Milestones B–D.
- **Categories dispatched:** `refactor_safety` (the step is a substitution and a directory move
  and must change nothing), `naming` (the step is a rename, and it turns on an owner's naming
  decision), `module_structure` (a module directory moved), `extras` (diff-matches-intent and
  report accuracy). `reliability` was **not** dispatched: A2 changes no test, adds no code path
  and moves no invariant — the coverage question A1's review answered has not changed, and the
  one forwarding test at risk was handed to `refactor_safety` as an explicit mutation mandate.

## 2. Verdict

**Approve-with-changes.** The substitution is *provably* pure — the `refactor_safety` agent
inverted the diff rather than reading it, and got the base file back byte-for-byte. Every
finding is documentation: two of them are statements the diff's own new prose made wrong or
incomplete, one is a leftover the step's verification could not see, and the rest are report
arithmetic and residue.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib` | 0 | 2,839 passed / 0 failed / 5 ignored — **unchanged from A1** |
| `cargo test --lib ng::` | 0 | 1,540 passed / 0 failed / 2 ignored — **unchanged from A1** |
| `cargo test --examples` | 0 | 52 passed / 0 failed |
| `cargo doc --no-deps` | 1 | 12 unresolved links, all pre-existing, **none** in the renamed module |
| `grep -rn "RecordReader\|record_reader" src` | 1 | no matches |

Four dumps byte-identical to the `8cf6f03` baseline by `cmp`; walk probe anchor exact.
Findings labeled "Needs verification": **0**.

**How byte-identity was established a second way.** The `refactor_safety` agent applied the
*reverse* substitution to all twelve changed files and diffed against `git show HEAD:<old
path>`: `read/filtering.rs` and `read/input/test_fixtures.rs` come back **byte-identical**,
`container.rs` never changed at all, the four reader files differ only by the added doc blocks,
and the remaining five differ only by rustfmt `use`/`mod` re-ordering the longer names force.
The only non-comment bytes in the whole diff that are not identifiers are two
`debug_struct("…")` labels, and nothing `{:?}`-formats those readers. On that evidence the
dumps' byte-identity is a consequence rather than a coincidence.

**The base occurrence counts were independently reproduced** by two agents: 11 / 9 / 17 / 30
bare / 23, per-file distribution matching file-for-file, with the new names appearing at base
only under `doc/` and never in `src/` — so nothing unrelated was captured. Of the 67 base
occurrences of the type names, none is a test name and none is an assertion string.

**The one forwarding test was re-mutated.** `the_enum_forwards_every_contract_method_to_its_arm`
still fails when the enum's `begin_region` body is replaced by `Ok(())` (marker `grep`-asserted
present before the run, per the `cargo fmt` hazard): `left: 0, right: 2`, on its own message.
The agent also mutated `other_sample_records` to a constant `0` — pinned too, though by
`a_shared_cram_serves_each_open_only_its_own_reads` rather than by the enum's own test.

## 4. Open questions and assumptions

1. **Should the four `pub(crate) mod` declarations inside `aligned_reads_reader/` become plain
   `mod`?** (affects Mi5) Nothing outside the directory names `bam` / `cram` / `in_memory` /
   `container`; the agent tightened all four and both `cargo check --lib` and
   `cargo test --lib --no-run` were clean. It is an internal visibility narrowing, not a rename
   — **deferred to Checkpoint A**, where it joins A1's Mi2 (the `pub` on
   `NoodlesRawAlignedRead`) as one visibility decision rather than two.
2. **How wide should the design-doc sweep be?** (affects Mi6) Five *live* ng design docs still
   name `record_reader/`, and A1's review found the same class. Deferred to Checkpoint A so it
   is done once for all three renames rather than three times partially.

## 5. Top 3 priorities

1. **M1** — the module's contract list states as a universal something the CRAM arm does not do,
   and A2's new prose is what makes that list read as complete.
2. **Mi1** — the module doc's inventory says the CRAM arm does not exist yet. It does, on the
   live path, and A2 rewrote that very sentence.
3. **Mi2** — the step's verification grep cannot see the old name spelled with a space; ten
   sites fell through, one of them a runtime panic message naming a type that no longer exists.

## 6. Findings

### Major

**M1: src/ng/read/input/aligned_reads_reader/mod.rs:69 — the contract list holds a universal with a live counterexample**
**Categories:** refactor_safety · **Confidence:** High

The list states: *"**Records come out raw.** In particular `read_group` is cleared, never
stamped."* `CramAlignedReadsReader::read_next` **does** stamp it (`cram.rs:189`,
`buf.read_group = Some(container.read_group(i));`), which the arm's own doc calls out as the
documented exception. The bullet is pre-existing text — but A2 added, fifty lines above it,
"**What every arm yields is undecoded** … so the module does", which is what makes the module
doc read as a complete statement of what an arm yields.

**Failure scenario:** this list is, by its own admission, the one place a new arm's author is
told to check their arm against, and it is prose that "does not fail a build". A universal with
a live exception is exactly the shape that produces an arm clearing a field it should have
carried — and a read group silently lost is a wrong per-read-group tally, not a crash.

**Fix:** amend the bullet, not the new paragraph — state the clearing rule, then the CRAM
exception and why it exists.

### Minor

**Mi1: src/ng/read/input/aligned_reads_reader/mod.rs:74-80 — the module doc's inventory is false**
**Categories:** module_structure, naming (convergent) · **Confidence:** High
"What is here so far: `InMemoryAlignedReadsReader` and `BamAlignedReadsReader`. The CRAM arm
lands in Milestone E" — but the enum has `Cram(CramAlignedReadsReader)`, `open_bam.rs:396`
constructs it on the live path, and `container.rs` is not mentioned at all. The diff rewrote
this sentence to substitute the new type names and left the claim wrong.

**Mi2: ten sites spell the old type name with a space, and the step's verification cannot see them**
**Categories:** naming, module_structure, refactor_safety (three-way convergent) · **Confidence:** High
The plan's check is `grep -rn "RecordReader\|record_reader" src`, which matches CamelCase and
snake_case only. `grep -rniE "record[ -]readers?" src/ng` finds ten live sites, including
**`bam.rs:148`'s `unreachable!` message** — *"a .crai index reached the BAM record reader"* —
which names a type that no longer exists and would be read by whoever hits that panic. Three
are in sentences this diff edited: `open_bam.rs:1799` now says "a record reader positions…
(`aligned_reads_reader/mod.rs`)", path renamed, subject not. **This matters beyond A2**: the
same grep would certify A3 the same incomplete way.

**Mi3: the enum's own naming justification contradicts the module's title and the impl report**
**Categories:** naming, extras (convergent) · **Confidence:** High
`mod.rs` argued *"The name says reads because that is what the file holds"* — while the module's
title line says "Where a cursor's **records** come from" and the impl report §6 defends keeping
"records" as "the layer's own vocabulary for **what a file holds**". That is the same premise,
claimed for both words. A2 also reworded the module doc's opener to "reads" while deliberately
leaving the contract list on "records", so the module carried three vocabularies in seventy
lines — and it was the *reword*, not the rename, that broke the agreement.

**Mi4: the impl report's arithmetic is wrong in two places**
**Categories:** extras, refactor_safety (convergent) · **Confidence:** High
(a) §3 says the undecoded statement landed "in **four** places" then enumerates five
(`grep -rn "yields is undecoded" src` → 5). (b) §3 says "the four type names, at all **67**
sites across eleven files: the enum, its three arms, **and the module path**" — 67 is the type
names only; the module path is a further 23, so 90 substitutions, and the 67 span nine files,
not eleven. §2's five counts are all correct and both agents reproduced every one.

**Mi5: the four `pub(crate) mod` declarations are wider than the surface needs**
**Categories:** module_structure · **Confidence:** High
Nothing outside the directory names `bam` / `cram` / `in_memory` / `container`. The agent
tightened all four to plain `mod`; `cargo check --lib` and `cargo test --lib --no-run` were
clean. `container` is the one that matters — a CRAM decode internal reachable crate-wide from a
directory whose doc says "one arm per format". **Open question 1.**

**Mi6: five live ng design docs still name `record_reader/`, and the renamed module cites back into them**
**Categories:** refactor_safety, naming, module_structure (three-way convergent) · **Confidence:** High
`arch/alignment_cursor.md` (the module tree at lines 27–32, plus ~14 more),
`spec/alignment_cursor.md:522`, `spec/alignment_file.md:40`, `arch/read_groups.md:24/286/343`,
`ng/README.md:86`. Not cosmetic: `aligned_reads_reader/mod.rs:7` and `bam.rs:30` cite back into
these, so the round trip lands a reader on a tree diagram naming a directory that is not there.
`read_groups.md:343` already carries one generation of this ("now `record_reader/bam.rs`") and
is now wrong a second time. Historical reports and the `PROJECT_STATUS.md` narrative are frozen
records and were correctly left alone. **Open question 2.**

**Mi7: two links in this plan's own spec and arch were broken by the `git mv`**
**Categories:** extras · **Confidence:** High
`spec/read_filtering_stages.md:147` — a live evidence pointer inside §4 "what will bite you" —
and `arch/read_filtering_stages.md:231` both point into `record_reader/`. Distinct from Mi6:
these are *this step's* governing documents, and the step is what broke them.

### Nits

Four "a `AlignedReadsReader`" where the substitution carried the old article through, while
`mod.rs` was hand-corrected to "An" — which is the tell that the others were not looked at; the
ASCII layer diagram in `region_records.rs` lost its column alignment (the new name is six
characters longer); two comment lines pushed past 100 columns, which `cargo fmt` will never
flag because it does not reflow comments; `cursor.rs:25`'s "the BAM arm lands at Milestone C and
CRAM at E", stale in the same way as Mi1 and on a line A2 edited; the two hand-written `Debug`
impls pick fields without destructuring `Self` (pre-existing, and `finish_non_exhaustive()` at
least declares the omission); the in-memory arm's new paragraph slips its pronoun between the
`Vec<RecordBuf>` handed in and the value yielded.

## 7. Out of scope observations

- `aligned_reads_reader/mod.rs`'s test doc says "Every method the enum offers reaches its arm"
  but exercises three of four. The `refactor_safety` agent **verified the fourth is pinned**
  anyway (constant-`0`-ing `other_sample_records` fails
  `a_shared_cram_serves_each_open_only_its_own_reads`), so this is a doc overclaim, not a test
  that cannot fail. Related: that test's `header()` assertion compares only
  `reference_sequences().len()`, the weak form `in_memory.rs:315-318` documents as having let a
  mutation through; the enum-level check was never strengthened the same way.
- `sample_cursor.rs:55` — `merged.cursors[0]` indexes with a literal on an invariant (`Merged`
  holds ≥2) that no type enforces. Pre-existing.

## 8. Missing tests to add now

**None.** A2 adds no code path, changes no invariant and moves no test. The one test whose
subject the rename could have hollowed out was re-mutated and still fails for its own reason;
the one method its doc overclaims about was separately shown to be pinned elsewhere. The
`reliability` category was deliberately not dispatched for this step, and nothing the other
four found suggests it should have been.

## 9. What's good

- **The reviewer inverted the diff instead of reading it.** Reverse-substituting all twelve
  files and diffing against the base is a complete proof for a rename step, and it is cheaper
  than the four dumps.
- **The CRAM arm's extra paragraph earns its length.** "Undecoded" genuinely has two senses
  there — `decode_container_at` really does build a `RecordBuf` per record before anything can
  be served — and the doc names both without muddling them. Independently verified by two
  agents against `container.rs:314`.
- **The substitution was ordered longest-name-first, and the report says why** — `RecordReader`
  is a substring of all three arm names, so the other order would have produced a name built
  from a half-renamed identifier.
- **`git mv` preserved rename detection** on all five files (similarity 100/95/94/93/86 %), so
  `git log --follow` still works.

## 10. Commands to re-verify

```
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --lib
cargo test --lib ng::
grep -rn  "RecordReader\|record_reader" src doc/devel/ng   # the widened check — Mi2, Mi6
grep -rniE "record[ -]readers?" src/ng                     # the blind spot — Mi2
```
Plus the four acceptance dumps and the walk probe against the `8cf6f03` baseline.
