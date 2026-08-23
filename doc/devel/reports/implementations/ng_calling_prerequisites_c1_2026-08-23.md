# ng — calling prerequisites, C1: the row that will carry a partial observation

**2026-08-23**, branch `ng-calling-prerequisites`. Step C1 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §5.4 and
[`arch/read_likelihoods.md`](../../ng/arch/read_likelihoods.md) §2.1.

**The type only, and the field it sits in is always empty.** A read that ran out inside a locus
says the sample carries *at least* what it saw; the merge throws that away today, and the calling
step has a censored term specified with nothing to read. This step gives the evidence a place to
be carried. The next one puts it there.

---

## 1. What changed and why

A partial observation is not a shorter allele. Its bases stop where its read's witness stopped, so
padding them out to the locus span and interning them would put an allele in the table that no
molecule carried — and it would read as a **short** allele, which is the one direction the model
must not be biased in ([`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §5.1). That
is why the merge's projection refuses one and why its iterator of projectable sequences exists.

What the calling design asks for is not to project them but to **keep them keyed by the stretch
they witnessed** (§5.4, corrected 2026-08-21), so a candidate can be scored against its own
projection restricted to those positions (§5.3). The evidence matters most where a locus is wide
on the reference: over half the overlapping reads are partial at a 60-base repeat tract, and an
allele longer than a read can only ever be witnessed partially (§5.4.1).

## 2. Changes made

**[`src/ng/run/cohort_merge/build.rs`](../../../../src/ng/run/cohort_merge/build.rs)**

- `PartialObservation` — the witnessed locus positions, the bases as the mint recorded them, the
  read group, a read count and a quality sum. One entry is one `(record, sequence, read group)`,
  **not one allele**: there is no allele.
- `SampleSupport` gains `partials: Vec<PartialObservation>`, filled with an empty vector. That is
  not a placeholder — it is what the merge does today, written down: both branches of the
  derivation drop every observation whose witness is not `Complete`.
- One test pinning the three things C2 must not move: a sample carrying a partial reaches the
  built locus with an empty row list, the partial's bases are not an allele of the table, and the
  locus's existence is unchanged. **C2 inverts the first and leaves the other two standing.**

## 3. Deviations from the plan

Two, both small, both recorded here rather than escalated.

**The plan says the row carries "the witnessed stretch (offset + length, off `ReadWitness`)".**
`ReadWitness::Partial` stopped being an offset and a length: it carries a canonical set of runs
(`WitnessedLocusPositions`), because the generic fold mints witnesses with holes in them and two
numbers can only describe a hole by swallowing it. The row carries the set.

**The runs are declared to be in the cohort locus's coordinates, not the record's.** The mint
measures a witness against the record it belongs to, and a cohort locus can hold several of one
sample's records and start before any of them. Doing the shift is C2's; saying which coordinates
the field is in is the type's, because a consumer indexing an allele's projection with an
unshifted run would read the wrong bases and nothing would say so.

**And one deviation from the architecture's sketch**, which types this as
`PartialObservation<'a>` inside the evidence view: the merge's row owns its bases and has no
lifetime, so the calling view can borrow `&[PartialObservation]` and there is one type rather than
two.

## 4. What the reviews changed

Three agents, each in its own worktree: what a locus with no partials must still see; every claim
re-measured; and test strength by mutation on a step that is deliberately behaviour-neutral.
**30 mutation runs between them.**

**The reversion test answers *yes* here, and that is the honest answer.** All three reverted the
change in whole and the module stayed green. For a type-only step that is correct rather than
alarming: what makes it non-vacuous is the compile error. `SampleSupport` has exactly one literal
construction, so C2 cannot forget the field.

**Two agents independently found the same unguarded rule, and it is the one this step's own
comment leans on.** The decision that a locus exists at all counts complete observations only —
the filter is `SampleLocusObservations::non_reference_and_compared_reads`, three lines in the
locus generator — and **it had no test at all**. Deleting the filter passed the merge's 246 tests
*and* the locus generator's 366, while turning a locus that nine partial reads leave quiet into
one that gets built. **This commit adds that test**, because it is the oracle the plan's own
verification table asks C2 for.

**A claim about what §5.4.2 settles was over-general, and it is the one a future reader would have
cited.** The first draft said locus existence is decided from complete observations "and stays
that way". §5.4.2 answers *no* under the heading "the merge's rule stays as it is **on the generic
path**", and then says the opposite for repeat tracts: a sample carrying an allele too long for a
read to span shows no complete observation, the filter reads that as *nothing varied here*, and
"one line of the rule has to change". This row is on both paths. Corrected, with the tract half
named as owed to whoever brings the STR path through the merge.

**A rule a consumer could act on was wrong.** The bases and the witnessed stretch were said to
have the same length "only where the read carried no indel over the stretch". They differ by the
*net* indel: a read carrying a two-base insertion and a two-base deletion inside the stretch comes
back with as many bases as positions and is still not a positional match. Anyone testing the
lengths for equality to decide "safe to index positionally" would be wrong on exactly that input.
Measured on a hand-built read, not reasoned.

**The test passed on a fixture with no partial in it.** Deleting the partial observation left all
its assertions true. It now asserts its own premise first. A separate check the reviews made —
turning the fixture's witness to `Complete` — makes the partial's bases a fourth allele, so the
allele assertion is load-bearing and is marked as such rather than looking like a duplicate of the
neighbouring test's.

**The ordering the field promises was under-specified.** "Ascending witnessed order" does not
order two entries that share a stretch, which is what a partially-witnessed substitution is. The
full key is now named — stretch, then read group, then bases — and C2 is told it owes the test the
sibling field already has.

**Three smaller corrections.** "A prefix or a suffix" predates the hole the generic fold now mints
and is corrected. "Never zero" was stated as a property of the type, which has public fields and
no constructor; it is the builder's to keep, as the sibling row's is. And what the row does *not*
carry is now recorded: four of the six quality sums a complete row holds, which strand bias, the
mapping-quality multi-mapper test and the read-position-bias term all read — at a repeat tract
those filters would see the complete reads only.

**What it costs, measured rather than waved at.** `SampleSupport` goes from **48 to 72 bytes**, a
50% widening, and an empty `Vec` allocates nothing (capacity 0, dangling pointer). Multiplied out
against the loci the merge holds at once — about one locus per 100 positions of resolved ground,
16 workers at the default region width — that is **768 bytes at one sample and 2.3 MB at 3,000**.
Against it: the observation cache the module's own docs call its dominant memory has a floor of
about **11.9 MB at 3,000 samples** before a single observation's heap. So this adds a real term,
under a fifth of a floor that is itself an underestimate, and no new scaling factor.

**One thing this step buys that the diff did not claim.** The module's determinism tests compare
whole `Debug` renderings, and its region-split test compares `SampleSupport` for equality. Both
were shown to read through the new field: two rows differing only in their partials compare
unequal and render differently. So the day C2 fills the field, those tests start covering it with
nothing new written.

**One recommendation not taken, and why.** Two reviews asked for a newtype over the witnessed
positions, because the field carries cohort-locus coordinates in the same type the mint uses for
record coordinates. The axis is in the field's name instead — `witnessed_in_locus` — and the
reason is recorded beside it: with one writer and no reader, a newtype would be minted for nobody.
It is the right move when a second consumer appears.


## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::run::cohort_merge` | **246 passed, 0 failed** |
| `cargo test --lib ng::locus_generation::tests` | **31 passed, 0 failed** |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 "redundant explicit link target" warnings — the same as the tree before |

**Amended after the commit: `cargo test --lib` is green in one run — 4,158 passed, 0 failed, 14
ignored, 597.13 s** — the 4,157 of the previous commit plus the one test this step adds.

The commit message says the whole-suite gate was split across runs, and at the time it was: two
runs on this tree stalled with the suite's two slowest tests outstanding — `joint::ssr_fit`'s
`a_drawn_stratum_returns_the_numbers_it_was_drawn_with` and `cohort_merge::serial`'s
`the_two_drivers_agree_on_random_layouts` — after the other 4,178 had passed, and each was then
verified green on its own. **The cause was the dev container, not the code**: the same suite ran in
621 s earlier in the session and degraded over hours of use. Stopping the container and letting
the next run start a fresh one restored it to 597 s and to a single clean pass. **Nothing about the
tree changed between the split evidence and the clean run.**

## 6. Follow-ups

- **Nothing reads `partials` yet.** The consumer is the evidence view of
  [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md).
