# Code Review: the census lives inside the psp — Milestone B

**Date:** 2026-09-09
**Reviewer:** rust-code-review skill (orchestrator), read-only category agents
**Scope:** the working-tree diff of each step of Milestone B, one review per step
**Milestone A's review:** [psp_census_pair_milestone_a_2026-09-09.md](psp_census_pair_milestone_a_2026-09-09.md)

---

## How this review was run

**The arrangement Milestone A settled on, unchanged.** The skill asks for one agent per category,
each in its own worktree; the build tree is 17 GB and the disk has 129 GB free, so that does not
fit. Instead one or two read-only agents cover grouped categories, forbidden to edit any file or
run `cargo`, and asked to **name** any mutation they want rather than run it. The orchestrator runs
the ones they name, serially, reverting each and re-running the module's tests on the restored
tree.

That trade paid three times in Milestone A: each step had at least one test that looked sound and
was green under the mutation it was written to catch.

---

## B1 — the verdict type

**Reviewed against:** the working tree over `b87c297d` — the new `src/ng/run/census_freshness.rs`
and two lines elsewhere. One read-only agent over five grouped categories: correctness and API
shape, the plan deviation, test strength, prose, and naming.

### Verification given to the reviewer

`cargo test --lib census_freshness` 4 passed. The tree's pre-existing red — 4 unformatted files,
4 examples that do not compile, 11 clippy errors in 6 files, 1 failing integration test — is
tabulated at the head of the
[implementation report](../implementations/ng_psp_census_pair_milestone_b_2026-09-09.md), and the
reviewer was told to ignore it unless the step added to it.

### Findings

**M1 — `Fresh` promised the whole judgement and the constructor gave it away for one comparison.**
`Fresh` is documented as a census in the format this build reads **and** recording the selection
this run rebuilds; `of_version` compared one `u16` and returned it. Spec §8 has
`regenerate-census` skip a psp on freshness *whole*, so a skip written against that doc would skip
psps it owes a rebuild.
*Fixed:* the constructor returns `Option<CensusVerdict>` and answers `None` — nothing to report —
for a matching version word, so no producer can reach `Fresh` through it. `Fresh`'s own doc now
says that how much it covers is the producer's to state, and that a judgement made from the psp's
head has looked at two of §4.2's three causes.

**M2 — `NoCensus` and `NotACensus` both described a trailer of zero bytes**, and nothing said
which wins, so B2 could implement either and both would read as correct against the docs.
*Fixed:* `NotACensus` is a trailer that *holds bytes* which are not a census, and the two causes
partition the psps.

**M3 — `AnotherFormat`'s two version fields could be written to disagree with the build that
printed them.** The variant is `pub` and so are its fields, so
`AnotherFormat { in_the_psp: 4, this_build_reads: 4 }` was one struct literal away — and `Display`
sent equality down the *not older* arm, rendering **"carries a census built by a newer version of
this program: it is version 4 and this build reads version 4"**. The suggested repair, private
fields, is not available: an enum variant's fields are public with the variant.
*Fixed by deleting the field.* The verdict carries the psp's version alone and the message reads
this build's from `VERSION`, which is the only place it exists.

**M4 — no test could see `Display` ignoring the field it was handed.** Both `AnotherFormat` tests
set `this_build_reads: 4`, which is `VERSION`, so a message printing the constant in both places
passed.
*Ran the mutation the reviewer named* — printing `{VERSION}` in the field's place — **and it was
green, 4 passed.** After the fix the same mutation fails two tests, because each expected string
names the found version as a literal and this build's through the constant.

**M5 — the doc counted three verdicts where there are four**, and then said "a command that only
has the first two available", which after the miscount does not say which two.
*Fixed:* four causes in three cost classes, and the cheap ones are named by what they read — the
footer, and the trailer's front — rather than numbered.

**Minor, all fixed.** The cost of reaching `AnotherSelection` was written as *seconds a sample*:
the reference read and the selection rebuild are once a run, and only a digest comparison is per
sample. *Reached the user as malformed* was past tense for something `decode_census` still does.
*The psp's records are whole* was asserted about a file whose blocks nothing had read — what the
open checked is the footer and the index. "Which is why they are one type rather than two" argued
the conclusion from the wrong premise, and is gone. `in_the_psp` / `this_build_reads` did not say
what they held; what survives is `version_in_the_psp`. `Hash`, `PartialOrd` and `Ord` are derived,
for spec §4.1's report grouped by cause, with a line saying that two samples at two different old
versions are two causes. And the module doc no longer uses the present tense for two commands that
are not built.

**`VERSION`'s doc contradicted its own value** — "2 since 2026-08-16" over `= 4`. Pre-existing, and
this step is what puts it on the public API surface. *Fixed* from the history rather than from
memory: 1 at creation on 2026-08-14, 2 on 2026-08-16, 3 and 4 both on 2026-09-05.

### Recorded, not fixed

**A census whose magic and version are this build's and whose sections are damaged has no
verdict.** No cheap read can tell it from a whole one; it reaches the user as the census reader's
own `CensusError` when something decodes it. The type doc says so rather than leaving the next
reader to find the gap.
