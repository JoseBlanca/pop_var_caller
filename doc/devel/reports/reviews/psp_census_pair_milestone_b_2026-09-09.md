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

---

## B2 — the cheap half of the judgement

**Reviewed against:** the working tree over `66427170`, three files — the psp reader, the census
file's own module, and the judgement. One read-only agent over five grouped categories:
correctness, reliability and errors, test strength, layering and API, and prose. **It named seven
mutations; six were run** (the seventh was superseded by a fix that changed what the code claims).

### Verification given to the reviewer

`cargo test --lib census_freshness` 9 passed, `--lib psp::reader` 59, `--lib census_file` 22. The
tree's pre-existing red is tabulated at the head of the
[implementation report](../implementations/ng_psp_census_pair_milestone_b_2026-09-09.md).

### Findings

**M1 — nothing in the suite failed if the judgement read the whole trailer**, which is the one
property the step exists to have. The reviewer worked it out by hand: every one of the five new
tests passes with `trailer()` in place of `trailer_head`, because the verdict is the same either
way and the verdict was all anything looked at.
*Ran it: 64 passed, no failure.*
*Fixed with the instrument the census's own reader already carries for the same argument.*
`census_file.rs`'s counting reader says it plainly — *an implementation that decoded a whole file
and handed back a slice would match every value a section-by-section reader gives and deliver none
of the memory the by-section design exists for; only the byte count tells them apart.* `PspReader`
now counts the trailer bytes it hands out, and two tests read the counter: ten bytes to judge a psp
carrying a census of a few thousand, and zero to judge one carrying none.

**M2 — the `# Errors` section named a failure that cannot reach this function.** It said a psp
*whose footer points past its own end*; `PspReader::open` proves the trailer ends exactly where the
footer begins (`reader.rs`'s footer check), so such a file is refused before it can be judged. The
sentence was load-bearing, because it is the argument for a read failure being an error rather than
a verdict.
*Fixed:* the two failures that can happen are named, and `open` is credited with ruling out the
third.

**M3 — the container module asserted what is in a trailer**, two lines below its own sentence
saying it must not. `psp_file_format.md` §3.4 keeps the payload opaque *so that* adding to it is a
writer-side change and not a container version bump, and calls that property worth more than any
list of payloads would be.
*Fixed:* the size argument stays — a payload can be tens of megabytes, ng's walk puts a census in
one — and the claim that the trailer *is* the census, and the census-shaped "first ten bytes", are
gone.

**M4 — `trailer_bytes as usize` wraps on a 32-bit target.** Pre-existing in `trailer()`, and this
step is where it now lives.
*Fixed:* `usize::try_from(…).unwrap_or(usize::MAX)`, a no-op on 64-bit.

**M5 — the shared fixture was reached through a `pub(crate) mod tests`.** The project's own answer
to *another module's tests need this fixture* is a `tests_support` module beside the tests —
`psp::writer` has one, six modules import from it, and this step's own tests already take
`a_header` from it. Exposing the whole tests module leaves a door open for anything later added
inside it.
*Fixed:* `tests_support` holds the census fixtures, `mod tests` glob-imports them and is private
again.

**M6 — "the psp's head" names the part of the file the function does not touch.** The format calls
the front of a psp its *header*; this function reads the footer and the front of the trailer, both
at the tail.
*Fixed:* `what_the_footer_and_the_trailers_head_say_about_a_census`.

**Minor, all fixed:** the module summary still said "the verdict and nothing that reaches it", which
this step falsifies; the function's own summary said *two short reads* where it performs one or
none, and left `NotACensus` out of the answers it lists; `trailer_head`'s parameter was `bytes`
where it is a count, and is `at_most`; `BYTES_THAT_NAME_THE_VERSION` sat 530 lines from the magic
and the version it is made of; the fixture helper's `sample` argument does no work until B4 and now
says so; and a version *newer* than this build's goes through the judgement and not only through
the message.

**Carried to B4, at the reviewer's prompting.** The judgement returns
`Result<CensusVerdict, PspReadError>`, and B4's contract is *every psp judged, no early return*. A
B4 written with `?` reproduces for read failures exactly the bug spec §4.1 exists to fix — one
unreadable psp in five and the user is told about one sample. B4 holds a result a sample rather
than propagating the first.

### The mutations, and what each showed

| mutation | outcome |
|---|---|
| the judgement reads the whole trailer | **green before the fix, 64 passed**; 1 test fails after |
| the head read ignores the trailer's length | 13 tests fail |
| the empty trailer is not answered out of the footer | 2 tests fail |
| the two damage verdicts are swapped | 4 tests fail |
| this build's own version word is not read as fresh | 2 tests fail |
| the version word is always read as this build's | before: 1 test failed, in the *other* module; after: 2, including `census_file`'s own |

**One the reviewer named stayed green, and the claim moved rather than the test.** Deleting the
empty-payload short-circuit in `trailer_head` fails nothing, because a zero-length read reads zero
bytes and the counter counts bytes. What that branch saves is a seek and a syscall; the doc now
claims *no byte of trailer is read*, which is measured, instead of *no read at all*, which is not.

### Confirmed by the reviewer and not changed

Every verdict is right for every input a psp can present, including a trailer of exactly the magic
with no version word behind it, and a version above this build's. `trailer()` behaves as it did
before this change for all four of its callers. And the 9-byte fixture's reasoning holds: the psp
format proves the footer begins where the trailer ends, so a reader that took ten bytes unclamped
would compose a version word out of one census byte and one footer byte.

---

## B3 — the cohort agrees on its settings

**Reviewed against:** the working tree over `89083b0b`, four files. One read-only agent over four
grouped categories: correctness and blast radius, error design, test strength, and prose. It named
six mutations; all were run, one of them in a rewritten form because the one named did not compile.

### Findings

**M1 — "nothing checked them until now", twice, and it is false.** `PspVariantCaller::open` has
compared every psp's catalog and criteria against the run's segmentation since psp mode was built,
so `call-from-psps` already refused a cohort typed two ways — naming one sample and the run rather
than the pair. What this step adds is the pair, and a refusal for the commands that never build a
run segmentation to compare against.
*Fixed in both places*, and the correction is what the implementation report leads B3 with, because
it changes what the step is worth.

**And the spec says the same wrong thing.** `psp_census_pair.md` §6's argument for the decision
opens with `first_difference` being "called only from its own tests"; it is called from
`psp_caller.rs`. The decision itself is unaffected. Left for the checkpoint rather than edited here.

**M2 — the new error's doc named commands that do not have the check.** `estimate-parameters` opens
its cohort through `open_census_cohort`, which compares the analysed regions alone, and
`regenerate-census` is not built. *Fixed:* it names `call-from-psps` and `generate-census`, and
points at plan step C2 for the third.

**M3 — the command-level test's account of its own fixture was wrong.** The second walk's purity
floor of 0.99 types no tract differently, because the fixture's reference is one base repeated and
every tract on it is perfectly pure. *Fixed:* the doc says what actually differs — the criteria
record in the header, which is what the opener compares.

**M4 — that test's assertions could not tell the new refusal from the old one.** Without the
opener's check the run reaches the caller's, which names the same field; only the sample pair
separates them, and the test asserted the two names separately.
*Fixed:* it asserts the sentence, *samples zeta and alpha do not agree on the set of repeat-tract
criteria*.

**Minor, all fixed:** a sentence that said a weak check "would pass this test" where it means the
test would fail, and pointed at a helper the test does not use; the module's opening summary still
said the opener settles "the ground they agree on"; the catalog test's doc argued an ordering the
test cannot see; the `MinCopies` fixture built a second segmentation to clone its inputs; the
command test re-walked both samples to obtain one psp; and `first_difference`'s own test repeated
the field-name literal that the routing constant exists to keep in one place.

**Minor, fixed with a test rather than prose:** the precedence between the two refusals — a cohort
differing in both the ground and the catalog is refused about the catalog — was inherited from
`first_difference` and pinned nowhere in the opener.

**Confirmed by the reviewer and not changed.** The three callers of `OpenPspCohort::open` are all
places a refusal is right. The 16 tests that depend on the `censuses_written_beside_the_psps`
fixture cannot build a disagreeing cohort: every one fills its directory from a single
`run_generate_psps` call, so the three settings are equal by construction. And
`PspVariantCaller::open`'s per-file loop still earns its place — the contig table it also checks is
not part of what the opener forces equal.

### The mutations

| mutation | outcome |
|---|---|
| the check reverted to the ground alone | 4 tests fail |
| the criteria dropped, the catalog kept | 2 fail |
| only the first psp examined | 5 fail |
| the catalog routed to the ground's refusal | 3 fail |
| the field-name constant's value changed | 1 fails, and only its own — which is what the constant is for |
| the caller's loop compares the first psp every time | **all 74 pass**, which settles that its segmentation check is now cohort-wide by construction |

The routing mutation as the reviewer wrote it — exchanging the two arms' bodies — does not compile,
since the variants carry different fields. A mutation that does not compile says nothing, so it was
run in the form above.
