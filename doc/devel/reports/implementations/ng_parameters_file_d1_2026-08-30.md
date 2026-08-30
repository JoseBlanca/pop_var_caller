# ng parameters file — D1: the bindings, derived from the run's own inputs

**Date:** 2026-08-30
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone D, step D1
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §3.1, §6
**Code:** [src/ng/calling/parameters_file/bindings.rs](../../../../src/ng/calling/parameters_file/bindings.rs) (new),
and the changed `of_run` in
[from_run_parameters.rs](../../../../src/ng/calling/parameters_file/from_run_parameters.rs)

---

## 1. What D1 turned out to be

The plan's D1 is one sentence — *the reference's content digest, the ordered sample list by name,
the read-group table, and the census recording terms the fit ran under — **written by B1 and read
by C2***. Three of those four were already written and read. The sample list and the read-group
table come straight off the run's own `ReadGroups`, and `validate` already refuses a file whose two
lists disagree with each other.

**What was missing is the step from the run's inputs to the two bindings that are not names.** A
reference digest and a census identity are both derived values, and until this commit neither was
derived anywhere: `of_run` took the digest as a `&str` — its own documentation said *"nothing here
can check that it is one"* — and took a `CensusIdentity` a caller had built by hand. Nothing in the
tree could build one. **D2 cannot compare a file's bindings against a run's until there is one
function that produces both sides**, which is what this commit is.

So D1 here is: **derive both bindings from the types the run already holds** — the census
identity minted here, the reference digest computed upstream and *spelled* here.

## 2. What was built

**`CensusIdentity::of(&RecordingTerms)`** — one term a value, in the words
`RecordingTerms::first_disagreement` uses and in the order it checks them: the seven selection
values first, then the loci actually kept, the per-stratum locus counts, the per-locus read cap,
the depth ladder edges and the per-position depth cap. Twelve terms.

Three choices inside it, each recorded because it could have gone the other way:

- **`RecordingTerms` is destructured without `..`**, which is that type's own convention upstream
  (`first_disagreement` and `SelectionTermsDigest::of` both do it). A thirteenth value added to the
  census stops this compiling rather than quietly dropping out of the identity — and a value that
  drops out is exactly the failure the binding exists to prevent.
- **The kept loci are digested whole *and* block by block.** The census compares both
  (`CensusLociDigest` derives `PartialEq` over `whole` and `per_block`), so an identity built from
  the whole alone would call two censuses the fit refuses to pool identical.
- **The two caps are digested like everything else**, though each is one small integer a person
  could read. Writing them as numbers would be friendlier, and that is the wrong direction for a
  binding: the file's own editing rule is that a readable key is an editable one, and a binding's
  whole use is an equality nobody should be able to satisfy by typing.

**A hex spelling, in one place.** `hex_digest` is 32 characters of lower-case hex, and it is the
only place either binding's text is produced — so the string a run writes and the string a later
run compares it against cannot be two spellings of the same bytes.

**`of_run` now takes the run's own `&ReferenceDigest`** and spells it here. The census still
arrives already minted, and **that asymmetry is deliberate**: a mismatched reference is refused
(spec §6), so the run's own digest is the only value `of_run` could ever write; a mismatched census
*demotes*, so a run that read a file fitted under other terms and writes its parameters out again
(spec §7) has to write back the terms it read. Which census a file names is its caller's to say;
which reference it names is not.

## 3. Two changes to the shared fixture, and why they are D1's

**The fixture's reference digest was eight bytes wide.** `A_REFERENCE_DIGEST` was
`"0123456789abcdef"` — 16 hex characters where a reference digest is an MD5 and spells as 32. Both
golden files carried it, so the artefact a reader learns the format from taught the wrong width. It
is now `THE_REFERENCE_A_RUN_FITTED_AGAINST`, a real `ReferenceDigest`, shared by every fixture in
the module.

**The fixture's census had one term.** No run can produce a one-term census identity, and the
golden file is the record of what the file looks like. It is now **a minted identity with its
digests replaced** — so the names and their order come from `CensusIdentity::of` and cannot drift
from it, and a thirteenth value added to the census fails an assertion rather than leaving the
fixture a term short.

**Only the digests are synthetic, and that is a choice against the alternative.** Pinning the
minted digests in a golden file would make tuning an unrelated default rewrite this module's
testdata: three of the twelve are taken over another module's `Debug` rendering of its own
defaults — `SelectionTermsDigest::of` hashes `format!("{catalog_built_under:?}")` and
`format!("{ssr_criteria:?}")`, and the depth ladder is digested from `DepthBinEdges::for_census()`.
The design-fidelity review checked that count and it is three.

## 4. One thing the writer gained

`[fitted_from.census]` had no note above it, where every other section of the file has one (step
B3). Twelve lines of hexadecimal with nothing above them is the section a reader is most likely to
stop at, so it now says what the terms are for, why a disagreement about one of them demotes
*every* number in the file, and what editing one buys. `reference_digest` gained a note too — what
the digest is taken over, which is what a reader needs to reproduce or check it. Both notes were
rewritten after the review; §9 says what they first got wrong.

## 5. What this does **not** do

- **It does not compare anything.** The three refusals are D2's and the demotion is D3's; nothing
  here reads a run's bindings against a file's.
- **It does not answer what a run with no census of its own should do.** Direct mode has no census
  (`run_streaming.md` §2 — no pre-pass, no psp), so there is nothing to compare the fourth binding
  against, and spec §6 does not say. That is F1's question, at the call site, and it is recorded in
  `PROJECT_STATUS.md` rather than guessed at here.

## 6. Tests

Six added, all in `bindings.rs`:

| test | what it holds |
|---|---|
| `every_term_is_named_as_the_census_names_it` | the twelve values `first_disagreement` checks, one edit each: the census reports a disagreement, the identity carries a term of exactly that name, and it is the **only** term whose digest moved |
| `every_part_of_a_composite_value_reaches_its_term` | the two values that are not scalars, part by part — the whole digest, a block's contig, its megabase, its own digest, an added block; a stratum's period, its repeats, its count, an added stratum. **Eight cases, added after review** (§9) |
| `the_terms_are_in_the_order_the_census_checks_them` | all 66 pairs: with two values moved, the census names the earlier one, and so does this identity. **The only test that can fail on an order** |
| `the_identity_names_each_value_once` | twelve terms, no name written twice |
| `a_digest_is_thirty_two_characters_of_lower_case_hex` | every minted term |
| `a_digest_spells_every_byte_as_two_lower_case_characters` | both ends of the byte range, in order — a dropped leading zero or upper case fails |

## 7. Validation

All in the container, by absolute path:

- `cargo fmt` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib ng::calling::parameters_file`: **145 passed, 0 failed** (139 before).
- `cargo test --lib`: **5,516 passed, 0 failed** (5,510 before) — the six new tests are the whole
  difference.
- `cargo doc --no-deps`: **25 unresolved-link errors**, the pre-existing baseline, and **23**
  `redundant explicit link target` warnings, also the baseline. *(An earlier draft of this line
  said 25 "unresolved intra-doc links"; the review measured that the unresolved ones are `error:`
  lines and the warnings are a different lint, and that the first draft of this step's own doc
  comment added a twenty-fourth of them. It uses a bare intra-doc link now.)*
- Both golden files regenerated by the two `#[ignore]`d helpers, and the diff read line by line:
  the reference digest widens from 16 to 32 characters in both, eleven census rows are added, and
  the written file gains the seven-line section note. Nothing else moved.

## 8. Mutation testing

**Fourteen mutants against `bindings.rs`, on the final tree; every one fails a test.** Each was
applied to a restored-pristine copy and the restore checked with `diff` before the next.
Command: `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`.

| mutant | tests it fails |
|---|---|
| drop the twelfth term | 121 |
| swap two terms' names | 7 |
| hex without the leading zero (`{byte:x}`) | 7 |
| hex in upper case (`{byte:02X}`) | 7 |
| rename one of the five hand-written terms | 6 |
| kept loci: the whole digest, no blocks | 3 |
| strata without their locus counts | 3 |
| kept loci: blocks without their megabase | 2 |
| kept loci: the blocks, no whole digest | 1 |
| kept loci: blocks without their contig | 1 |
| kept loci: blocks without their own digest | 1 |
| strata without their period | 1 |
| strata without their reference repeats | 1 |
| the strata as one grand total | 1 |

**Seven of those fourteen survived the first draft of this step** — every one of them a *part* of
one of the two values that are not scalars. See §9.

## 9. What the review found

Three agents in isolated worktrees: correctness, design fidelity, and a geneticist reading the
produced file and forbidden to open a source file.

### The Major: the two composite values were tested for existence, not for contents

Ten of the census's twelve values are scalars and the twelve-edit table moved each of them
exactly. The other two are not — the kept loci are a whole digest plus a list of blocks, and the
per-stratum counts are a list of keys and counts — and each got **one** edit that moved **one**
part. So dropping the whole digest, a block's contig, a block's own digest, a stratum's period,
its repeats, its count, or replacing the whole count table with its grand total were all
invisible: seven mutants, and the suite said `144 passed` in every one of those worlds.

**What that would cost is the thing this step exists to prevent.** With the contig undigested,
two runs whose kept loci differ only in which chromosome a megabase sits on — an everyday
difference — mint byte-identical census identities. The fit calls those two censuses unpoolable;
the file would call them the same, and D3 would not demote.

Fixed by `every_part_of_a_composite_value_reaches_its_term`, eight cases, and by two fixture
changes without which the cases prove nothing: the counts start non-empty (an empty table makes
*any* surviving byte move the digest), and the blocks are two on different contigs.

**One of the eight cases was wrong on its first run and had to be measured to see it.** The
period edit moved a stratum's period from 2 to 9, which sends it past its neighbour in
`iter_sorted`'s order — so the bytes moved whether the period was digested or not, and the
"strata without their period" mutant still passed. At period 1 the sorted order is unchanged and
only the period differs. **This is why the mutant was run rather than the test being read.**

### The order was asserted and unpinned

`CensusIdentity::of`'s doc claimed the twelve are in `first_disagreement`'s checking order, and
no test compared the two: moving one value at a time cannot see an order, so any permutation
passed. The order is load-bearing — `first_disagreement` reports the *first* value two censuses
differ on, so when two have drifted the order decides which term a run names. Now
`the_terms_are_in_the_order_the_census_checks_them`, over all 66 pairs.

### Four doc comments that said something untrue

- **The module header claimed the digesting was "the census's own choice one level down".** It is
  for seven of the twelve; the census *file* writes the per-stratum locus counts, the read cap
  and the depth cap as values (`census_file.rs`, `encode_header`). The header now says which are
  inherited and which are this file's own choice.
- **It claimed four bindings.** Two are minted or spelled here; the sample list and the
  read-group table are names read off the run's own table.
- **It said both were "minted here".** The census identity is; the reference digest is *spelled*
  here and computed upstream by `ReferenceDigest::of`.
- **A fixture claimed "every value differs from every other"** where two of the seven selection
  values were both `StrRepeatCriteria::default()` — so a digest taken over the catalog's criteria
  where it meant this run's would have been invisible. They are different values now.

### What the geneticist could not answer

- **`reference_digest`: no algorithm and no subject.** Thirty-two hex characters could be MD5 or
  half a SHA-256, over the whole file or the bases or the contig list — and that reader keeps a
  soft-masked, an unmasked and a re-wrapped copy of the same tomato assembly. A refusal naming
  two strings would leave them unable to tell which of the three the file was fitted on. The key
  now carries a note: the MD5 of every contig's bases, uppercased, run together in file order —
  so soft-masking and line width do not change it and contig order does.
- **"Nothing here can be usefully edited" was the wrong claim**, and they said why: there is one
  very useful edit, which is to type the run's own digests in and have a file fitted on somebody
  else's evidence reported as fitted on yours. The note now says what an edit buys rather than
  that it is futile, and extends the warning to `reference_digest`, `samples` and `read_groups`
  above, which had none — and where the failure is the worse one, since those refuse.
- **Why *every* number is demoted for a disagreement about one line.** Stated, not explained. The
  note now gives the reason: the numbers were fitted together out of one store of evidence.
- **⚑ Where the run says which line differed — still unanswered**, and recorded in
  `PROJECT_STATUS.md`. The note promises it, D3 will produce it, and what an output prints from
  it is the emission step's document, unwritten (spec §11).

### Recorded and not fixed

Seven findings reach outside this step and are in `PROJECT_STATUS.md` for the owner: spec §12's
third question is closed by what shipped; a `.fai`-only run cannot produce the `ReferenceDigest`
this signature now demands, against §7's unconditional write; the census's own `reference digest`
term can never be the line a demotion names; the two caps are digested where the census file
writes them as numbers; and four comments in sections this step did not touch say things that are
false or backwards.

### One thing about the review harness itself

**`./scripts/dev.sh cargo test` run from inside a worktree tests the main repo, not the
worktree.** `scripts/dev.sh` sets `PROJECT_DIR` from its own location and passes `-w
"$PROJECT_DIR"`, so the container's working directory is the main checkout whatever the host cwd
was — a worktree holding a mutation reports the unmutated result. Two of the three agents hit
this; both noticed, because a run with no `Compiling` line finished in 0.32s. The form that works
is `dev.sh bash -c 'cd <worktree> && CARGO_TARGET_DIR=<worktree>/target-container cargo …'`.
