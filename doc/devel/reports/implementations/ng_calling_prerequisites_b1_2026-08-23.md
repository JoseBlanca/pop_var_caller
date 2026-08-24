# ng — calling prerequisites, B1: the merge's per-allele support gains a read-group axis

**2026-08-23**, branch `ng-calling-prerequisites`. Step B1 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §2.3 and
[`arch/read_likelihoods.md`](../../ng/arch/read_likelihoods.md) §2.1.

**One row per `(allele, read group)` where there was one row per allele.** 235 of the merge's 236
existing tests pass untouched, and a one-read-group sample's rows come out byte-identical —
measured on four locus shapes including composed reads and unified placements, not inferred. The
236th is the only pre-existing fixture with two read groups in it, and this step inverts what it
asserted.

---

## 1. What changed and why

When the merge records what a sample's reads lend an allele, it summed every read that reached
those bases. A read likelihood may fold an observation's reads into one term only if every one of
them would get the same number, and two reads showing the same bases from two lanes have
different error rates — so that sum crossed a boundary it may not cross
(`spec/read_likelihoods.md` §2.3). The locus generator already keeps the groups apart:
`SequenceObservation` carries `read_group` as part of an observation's identity. The merge was
throwing it away, and its own doc comment had booked the change as owed.

## 2. Changes made

**[`src/ng/run/cohort_merge/build.rs`](../../../../src/ng/run/cohort_merge/build.rs)**

- `SupportedAllele` gains `read_group: ReadGroupId`. `SampleSupport::supported` is now one row
  per `(allele, read group)`, in ascending `(allele, read group)` order.
- The per-sample tally was a `Vec<AlleleSupportTally>` indexed by allele. A pair cannot index a
  `Vec`, so it becomes a `Vec<AlleleGroupTally>` carrying its key, found by linear scan and
  sorted once per sample. **The scratch buffer keeps its allocation, which is what it was for**;
  what it loses is the addressing, and at an ordinary locus the scan looks at one or two entries.
- `AlleleBacking::read_group()` — which group a piece of evidence came from, in one place. A
  read has one read group because the group belongs to the library the fragment was prepared in,
  so a read seen at several of the sample's records has one too, and the first sighting is the
  answer. **B2 owns the assertion that the sightings agree**; this step reads the first because
  there is nothing else it could be.
- `SampleSupport::support_for` → `pooled_support_for`, now folding every matching row rather than
  finding the first, via a new `AlleleSupport::added_to`. **The name carries the warning**: a read
  likelihood must not use it, for the reason above. All 17 call sites are tests.
- Two doc comments that booked this change as owed are rewritten. The old text gave the reason as
  the STR path's — stutter is fitted per read group — and that reason still holds; what changed
  is that the generic path turned out to need the split too, so it is no longer owed to a later
  step.

## 3. The split between B1 and B2, which is not quite the plan's

The plan gives B1 "the type change" and B2 "attribution stops at the boundary". Taken literally
that leaves an intermediate commit in which the row carries a read group but the tally still
pools, so a two-lane sample gets one row labelled with one of its two groups — a knowingly wrong
output, committed on purpose. **So B1 here is the type *and* the key**: the rows it emits are
correct for every sample whose alleles come from single-record evidence, which is the ordinary
locus.

**B2 keeps the half that is genuinely subtle** and that the plan's own note is about: the
divided-read path, where one read is seen at several of the sample's records and its share of
each is apportioned. That is where a read's group has to be established and asserted, and where a
wrong answer is a quietly wrong `q_sum` rather than a crash. Recorded here rather than escalated,
as a deviation that keeps the plan's intent and its commit boundary.

## 4. What the reviews changed

Three agents, each in its own worktree: what a one-read-group sample must still see; every claim
re-measured; and the strength of the tests under mutation. **63 mutation runs between them.**

**The failure shape the two preceding steps had is absent here**: the change cannot be reverted
with the tests green. Dropping the read group from the tally's key reddens both new tests, and
dropping the sort reddens them too, each by a different assertion.

**The one severe finding is closed, and I checked rather than assumed it.** Reading a composed
read's group off the record's *first* observation instead of the one the read was sighted at
survived every test in the library — it files a read under a lane that produced none of its reads,
which is exactly what §2.3 forbids. The divided-read fixture added for the rounding grain turns
out to kill it: read 11 collapses into lane 1 and the row count drops from two to one. Applied the
mutation and ran it to be sure.

**The parity claim was weaker than it read, and the report said it plainly wrong.** "236 tests
pass unchanged" invited the reading that 236 one-group fixtures folded to today's shape. 235 did.
The 236th, `support_is_summed_within_an_allele_and_never_across_them`, is a two-lane fixture whose
doc said the two groups' reads "are summed" and whose assertion message said so too — and it
stayed green only because the renamed accessor re-adds across groups. Renamed, corrected, and
given row-level assertions so its green now rests on the merge rather than on the accessor.

**Two behaviours the tests did not defend.** Deleting the sort left 236 of the 238 green: the
tally used to be a vector indexed by allele, so ascending allele order was free, and it is now a
line a careless edit can remove. A new one-group fixture emits its alternative before its
reference and so would come back out of order without it. And `AlleleBacking::read_group`'s
multi-record arm had no test at all — both new tests built a single record — so reading the group
off the wrong sighting changed nothing anyone could see.

**A grain changed and neither the code nor this report had noticed.** Finishing a row rounds a
divided read's four count-like sums, and that happened once per allele; it now happens once per
`(allele, read group)`. So two lanes each holding half a forward read now round to one each, and
the rows total two where the single row would have held one. Neither answer is wrong — a row has
to be whole reads to be usable on its own — but the first draft of this step asserted the opposite
in a test *name*, `…is_two_rows_that_add_back_to_one`, on a fixture that could not see it because
its sample had one record. The test is renamed, a second test pins the real grain at one forward
read, and `AlleleSupport`'s doc says which grain it now rounds at.

**Three prose corrections.** The inherited byte figure — "entries of 32 bytes" — was the size of
the support and not of a row; a row was 40 bytes and is now 48, and the sparse shape's saving
against a dense record is gone at the very locus the same paragraph names (96 bytes either way at
2 of 3 alleles). What still pays is the cohort-wide case, and the paragraph now says that instead.
The reason given for a scan rather than a map was allocation, which a map cleared per sample would
not do either. And "most samples of most runs" had a measurement it was not quoting: 157 of 1,707
samples in a surveyed tomato archive carry more than one library.

**Six design documents stated the split as still owed** and are corrected in this commit:
`arch/read_likelihoods.md` §2.1 and its reconciliation table, `spec/read_likelihoods.md` §2.3 and
its §7 ownership table, and `spec/cohort_merge.md` §4.2 in two places. The merge's own comment had
booked the change as the STR path's to make; the generic path needed it first.

## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::run::cohort_merge` | **240 passed, 0 failed** — the merge's own suite |
| `cargo test --lib` | **4,151 passed, 0 failed, 14 ignored**, 576.61 s |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 "redundant explicit link target" warnings |

The tree before this step gave **4,147** and the same 24 and 12, so the count moved by exactly the
four tests added and this step introduces no doc-link failure. **It did introduce one and it was
caught here:** a doc comment linked to a `#[cfg(test)]` item, which does not exist in a doc build,
taking the count to 25. Named rather than linked now.

The merge's own suite held 236 tests before this step and holds 240 after, with 235 of the 236
untouched.

## 6. Follow-ups

- **B2** — the assertion that one read's sightings agree on their read group. Until it lands,
  reading the group off a read's *last* sighting rather than its first passes every test, and only
  a fixture where one read's sightings name different groups could tell them apart — which is the
  state B2 exists to rule out.
- **Nothing outside the merge consumes `supported` yet**, so the axis costs no other module a
  change today. The evidence view that will consume it is
  [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md)'s, and the arch
  already describes it as one entry per `(allele, read group)`.

- **The scan is unmeasured where it could stop being free.** The callback fires once per read on
  the multi-record branch, so a deep sample whose reads compose many distinct alleles does about
  45,000 comparisons at a locus where the old index did 300 lookups. At the one or two entries an
  ordinary locus holds it is nothing. Recorded beside the code; the fix if it ever matters is a
  sorted insert rather than a map.

- **Determinism has never been checked with more than one read group.** The module's four
  determinism tests compare the whole `Debug` rendering of a region's outcome, so they do now
  cover `read_group` — but every fixture behind them uses group 0, because the shared `sequence`
  helper hardcodes it. It holds by the uniqueness of the sort key rather than by measurement.

- **Two impl-plan documents record the merge's own steps without the axis** —
  `impl_plan/cohort_merge.md` step B3 and `impl_plan/calling_read_likelihoods.md`'s dependency
  note. They are dated records of completed and future steps rather than live design, so they are
  left as written; the specs and architecture they rest on are corrected.
