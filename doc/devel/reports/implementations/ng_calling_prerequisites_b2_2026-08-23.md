# ng — calling prerequisites, B2: a read's sightings must agree on their read group

**2026-08-23**, branch `ng-calling-prerequisites`. Step B2 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §2.3 and
[`arch/cohort_merge.md`](../../ng/arch/cohort_merge.md) §5.

**The merge now refuses a read that two of one sample's records name in two different read
groups, instead of filing it under whichever group its first sighting happened to carry.** Three
tests, one check, one rewritten doc comment; the merge's own suite goes from 240 to 243 and
nothing else moves.

---

## 1. What changed and why

B1 gave the merge's per-allele support a read-group axis: one row per `(allele, read group)`. A
read composed across several of a sample's records has to be filed under one of them, and B1 read
the group off the read's **first** sighting, on the argument that a read has one read group
wherever it is sighted — the group belongs to the library the fragment was prepared in. Nothing
checked it, and B1 recorded the check as owed in its own doc comment ("B2 asserts that they
agree") and in its follow-ups.

**A disagreement is one of two upstream defects, and both are silent.** Either the mint gave two
reads one chain id, so what the merge composes across records is not one read and its bases belong
to no fragment; or a read's group was stamped from the record rather than from the read. Neither
crashes. What comes out is a row whose `q_sum` was summed under a lane that produced none of those
reads, which is the pooling `spec/read_likelihoods.md` §2.3 forbids, at a locus that looks
ordinary.

The mint stamps the group from the read
([`open_record.rs:2093`](../../../../src/ng/locus_generation/pileup/open_record.rs) —
`read_group: active.read.read_group`), so the invariant holds today by construction. This step
makes it hold by check.

## 2. Changes made

**[`src/ng/run/cohort_merge/build.rs`](../../../../src/ng/run/cohort_merge/build.rs)**

- `AlleleBacking::read_group`'s multi-record arm compares every sighting after the first against
  the first and refuses on a disagreement, naming the read and both record regions. Release-level,
  like the chain-id check in `alleles_of_sample` a few lines above, and booked the same way: **when
  observations are decoded from a psp file this becomes corrupt input and must become a
  `RunError`** (`arch/cohort_merge.md` §5).
- The doc comment says what a violation would mean and why there is no answer to fall back to,
  replacing the sentence that booked the check as owed.
- Three tests. A three-record fixture in which one read's sighting at a chosen record carries a
  second read group: one test puts the disagreement at the **middle** record and one at the
  **last**, because a check that looked at only one sighting after the first would pass whichever
  of the two it happened to look at. The third builds the same fixture with one lane throughout
  and asserts it yields one row — so the two refusals are evidence about the disagreement and not
  about the fixture.

## 3. Deviations from the plan

None beyond B1's, which stands: B1 took the type *and* the key, leaving B2 the check on the
divided-read path. That split is recorded in
[B1's report](ng_calling_prerequisites_b1_2026-08-23.md) §3.

## 4. What the reviews changed

Four agents, each in its own worktree: what a legitimate sample must still see; every claim
re-measured; test strength by mutation; and whether the check covers the invariant it names and
sits in the right place. **45 mutation runs between them.**

**The change cannot be reverted with the tests green.** All four tried it: deleting the loop and
restoring the first-sighting read fails 2 of the module's tests, both through
`CohortObservation::over` rather than a helper. The failure shape that appeared twice earlier on
this plan does not recur.

**The reason given for the check was wrong, and the correction is the substance of this step.**
The first draft said a disagreement means one of two upstream defects: two reads given one chain
id, or a group stamped from the record. The second **cannot happen** — a record has no read
group; the mint stamps it from the read (`open_record.rs`, `read_group: active.read.read_group`).
And the reachable mechanism is neither: **a fragment's two mates are collapsed onto one chain id
on their name alone** (`chain_id_allocator.rs`, `pending_mates`), while the read group is
resolved from each SAM record's own `RG` tag (`read_groups.rs`), and nothing requires two mates
to carry the same one. So a well-formed BAM whose mates disagree about their `@RG` — or a merged
file in which two libraries reuse a read name — reaches this assertion with nothing defective
anywhere. Two agents found this independently, each tracing the chain from the tag to the panic.

**That makes the failure mode wrong by the architecture's own criterion, and it is recorded
rather than repaired here.** `arch/cohort_merge.md` §5 puts a fact about the data on the counted
side and keeps panics for bugs in whoever hands the work out; the release profile aborts, so a
header the user wrote can end a run. The check is kept as a panic because this module has no
`RunError` to raise, and §5 now carries both this variant and the chain-id one as owed at the
psp-decode step, together with the better repair: carry the read group on the pending mate and
refuse a mismatched second mate where the chain id is handed out. **That repair also covers a
form this check cannot see** — where such a fragment's mates *overlap* and agree, the walk sums
their base quality and gives it to one of them (`resolve_mate_overlap_at_pos`), stamping one
library's quality with the other's group in a single observation, with no second sighting to
compare.

**The message could not be acted on and is rewritten.** It named a chain id, which is comparable
within one sample only, without naming the sample — every neighbouring assertion in the file
names it — so `read_group` now takes the sample index. And it listed the two record regions in
the opposite order to the two groups printed beneath them, so pairing them positionally got both
backwards; each group is now named beside its own record, in that order.

**Three test gaps, each proved by a surviving mutation, each closed.**

- **Two records — the smallest and commonest multi-record sample — had no test.** Skipping the
  check for two-sighting reads passed all tests and returned a row with the whole read filed
  under one lane. Closed by a two-record fixture.
- **The baseline sighting itself was unpinned.** Taking the group off the *last* sighting instead
  of the first passed everything: at a disagreement after the first record the check still fires,
  and at a disagreement *at* the first record it quietly files the read under whatever the others
  carry. Closed by a fixture whose odd group is at the first record.
- **Which record the message blames was unpinned.** Swapping the two regions changed the output
  and no test noticed. Closed by making each `#[should_panic]` expect the whole clause, regions
  included.

**One survivor is left standing knowingly.** Downgrading the assertion to a `debug_assert!` passes
every test under `cargo test --lib`, because the test profile arms debug assertions and the
release profile does not. It is not invisible: `cargo test --release --lib ng::run::cohort_merge`
fails the refusal tests under the downgrade and passes them as written, which is how the
"release-level" claim was checked here. The gate as specified does not run that, and the same
hole covers the chain-id assertion this one is modelled on.

**Two claims about cost were wrong and are corrected.** "Two or three records at the loci that
reach this branch" had no measurement behind it and the module's own fixture contradicts it: the
generic mint writes a record at every covered position, so a sample holds as many records as the
locus is wide — six inside a six-base locus in `serial.rs`'s minted fixture — bounded by
`MaxCohortLocusSpan`, 50 by default. And a wrong group does not spoil one `q_sum`: the group is
half the tally's key, so the read's whole share moves. **What the cost actually is** was measured
by the fourth review at the deep end of the committed range — 50 records, 300 reads, 14,700
comparisons per sample per locus — at **+3.7 µs on a 320 µs locus build, against a run-to-run
spread of ±13 µs**: not separable from noise, because the same locus already sorts every sighting
and composes and divides every read across every record.

**One thing the check does not cover, now said in the doc comment.** On the single-record arm the
same invariant holds by construction — the read group is part of the key the mint buckets on —
but by the time the merge sees the observation there is one group and a list of chain ids, and
the per-read groups are gone. A psp file corrupted that way is undetectable there.

**A B1 follow-up is closed rather than merely untested.** B1 recorded that reading the group off
a read's last sighting rather than its first passed every test. After this step the two anchors
are equivalent by construction: either they agree and the answer is the same, or the check fires.

## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::run::cohort_merge` | **245 passed, 0 failed** — the merge's own suite |
| `cargo test --lib` | **4,156 passed, 0 failed, 14 ignored**, 644.37 s |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 "redundant explicit link target" warnings |

The tree before this step gave **4,151** and the same 24 and 12, so the count moved by exactly the
five tests added and this step introduces no doc-link failure. The merge's own suite held 240 and
holds 245.

**One gate beyond the stated set was run**, because the stated set cannot see the property the doc
comment claims: `cargo test --release --lib ng::run::cohort_merge` compiles with
`debug-assertions = false` and so is the only run that distinguishes this check from a
`debug_assert!`. It reports **245 passed, 0 failed** — the refusals fire with debug
assertions off.

## 3b. Documents corrected in this commit

- **`arch/cohort_merge.md` §5** gains the two assertions owed to `RunError` at the psp-decode
  step, and says why the read-group one is the sharper case — it can be reached by a well-formed
  BAM.

## 6. Follow-ups

- **The real repair is upstream and is not this check.** `PendingMate`
  (`locus_generation/pileup/chain_id_allocator.rs`) should carry the read group, and the
  second-mate branch should refuse a mismatch as a `WalkerError` — one comparison per paired read,
  on an error channel that already exists, and it covers both the case this check catches and the
  overlapping-mates case it cannot see. Out of this step's scope: it is the locus generator's,
  and it is a design decision about a failure mode rather than a coding slip.
- **The gate set cannot see a `debug_assert!` downgrade.** Adding
  `cargo test --release --lib <module>` for modules carrying release-level assertions would close
  it, for this check and for the chain-id one.
- **The panic's blast radius under parallelism is what it already was for the assertions beside
  it**: rayon joins every builder in the round before re-raising, and the cache is left advanced
  and not rewound (`parallel.rs`,
  `a_builder_that_panics_leaves_the_cache_in_the_callers_hands_and_advanced`).
