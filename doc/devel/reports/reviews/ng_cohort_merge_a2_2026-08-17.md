# Review — ng cohort merge, step A2 (the two derivations)

*2026-08-17, branch `ng-cohort-merge`, working tree at stash commit
`2e0f8adf`. Four category checklists, three sub-agents, two of them mutation-testing in
isolated worktrees. Per-category audit trail:
`tmp/review_2026-08-17_ng-cohort-merge-a2/`.*

## 1. Scope

- **Reviewed:** `SequenceObservation::matches_reference`,
  `SampleLocusObservations::reach`, `::non_reference_reads`, their tests, and the moved
  call site in `CensusWriter::add_generic`.
- **In-scope files:** [`src/ng/locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs),
  [`src/ng/parameter_estimation/joint/census.rs`](../../../../src/ng/parameter_estimation/joint/census.rs).
- **Categories dispatched:** `reliability` and `refactor_safety` (the two that can be
  answered by mutation, and the step's risk is a change to shipped code), `naming`,
  `smells`. **Not dispatched:** `defaults` (no parameter), `unsafe_concurrency` (none),
  `errors` (no error path), `tooling`, `module_structure` (no tree change), `idiomatic`
  (folded into `smells`).

## 2. Verdict

**Approve-with-changes.** No Blocker. Four Major, every one a real gap; all applied.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::locus_generation::tests` | 26 passed, 0 failed (at review time) |
| `cargo test --lib ng::parameter_estimation::joint` | 103 passed, 0 failed, 501.66s |
| `cargo clippy --all-targets --all-features -- -D warnings` | fails, 49 pre-existing errors, none in scope |

**Mutation testing, both worktree agents' numbers as reported:**

| agent | mutants | survived | changed-no-behaviour |
|---|---|---|---|
| reliability | 12 distinct, 13 runs | 3 | 0 |
| refactor_safety | 5 | 1 | 0 |

The two agents' survivor sets overlap on the containment mutant, found independently.
Every survivor was proven to change behaviour on a constructed fixture before being
recorded — neither agent reported an equivalent mutant as a finding.

## 4. Findings

### Major

**M1 — the predicate has no guard against bases of a different length from the
reference, and that is the whole indel half of what this caller does.**
**Categories:** reliability, refactor_safety — convergent, found independently.
**Confidence:** High, measured.
Replacing `*self.bases == *reference_bases` with `reference_bases.starts_with(&self.bases)`
left **58 census tests and 360 `locus_generation` tests green**; so did
`self.bases.starts_with(reference_bases)`, and so did `ends_with`. The mutants are not
equivalent: against reference `ACGT`, a deletion carrying `AC` is non-reference under
equality and reference under the first; an insertion carrying `ACGTT` is non-reference
under equality and reference under the second. Both agents proved the divergence with a
probe before recording the survivor. No fixture anywhere — in either caller's tests —
compared bases of a different length against the reference.
**Why it matters:** `non_reference_reads` is what the merge's keep rule sums. Under
either wrong implementation an indel stops counting, indel-only loci fall below
`min_alt_obs`, and nothing downstream can recover them: wrong results, no panic.
**Fixed:** three assertions — a trailing deletion, and insertions on both sides —
added to `matches_reference_compares_the_bases_it_is_given`, plus the reason on the
test. **Verified after the fix:** each of the two containment mutants now fails exactly
that test (28 passed, 1 failed).

**M2 — two more open-coded copies of the predicate ship in library code, while the doc
declared itself the only one.**
**Categories:** refactor_safety, smells — convergent. **Confidence:** High.
[`depth_and_alt_reads.rs:82`](../../../../src/ng/parameter_estimation/generic/depth_and_alt_reads.rs)
and `:234` both read `if *observation.bases != *locus.reference_bases`, over
`complete_observations()`, summing `num_obs` — not merely the predicate but the whole of
`non_reference_reads`, reimplemented twice, and live (`accumulators.rs` calls all three
entry points). Mutating the shared predicate to `true` left them untouched, which is how
the agent proved they never reach it. Two further spellings sit in
`locus_generation/pileup/tests.rs`.
**And the collision runs both ways:** that module's own doc claims *"the only place that
decides what counts as an alternative read, which is why it is its own file rather than a
method on the locus type"*, from `arch/parameter_prepass_generic.md` §2.3. One document
says the definition must not live on the observation type; A2 put it there and declared
itself unique. A reader arriving from either side is told the other does not exist.
**Fixed, by separating the two claims rather than choosing between them.** All four
remaining spellings now call the predicate; `matches_reference`'s doc claims only to be
*the one place the comparison is written*, and says explicitly that which observations to
ask about — the subset, the depth cap, the read-group grain — is still the pre-pass's.
Both documents are then true. **The arch wording is still worth one line from the owner**
(raised at Checkpoint A).

**M3 — the recorded cost of the complete-only rule had no test, in the direction the
existing test cannot see.**
**Categories:** reliability. **Confidence:** High.
`a_partial_that_agreed_with_the_reference_is_not_counted_against_it` covers a partial
that *agreed*, where the answer is 2 whether partials are excluded or compared properly
against their own stretch — so it separates the naive bug from today's rule, but not
today's rule from a correct partial-aware one. The behaviour that will surprise a reader
— a partial that *disagreed*, with 7 reads, contributing 0 — was untested.
**Fixed:** `a_variant_seen_only_by_partial_reads_is_not_counted`, which is what will fail
when step 7's censored likelihood changes the decision.

**M4 — `reach`'s doc claimed agreement with production "everywhere", and that is false at
the one input its own test exercises.**
**Categories:** reliability. **Confidence:** High, measured.
Production saturates the addition before subtracting, so a one-base region at `u64::MAX`
reaches `18446744073709551614` under its expression and `18446744073709551615` under
this one. **This code is the right one** — the last base of that region is `u64::MAX` —
so the defect was in the sentence.
**Fixed:** the claim is narrowed to "every region below the top of the coordinate space",
the divergence is stated with its direction, and
`a_locus_at_the_coordinate_ceiling_reaches_its_own_end` now asserts
`production_reach(u64::MAX, 1) == u64::MAX - 1`, so nobody restores "agreement" by
reintroducing arithmetic that both panics and answers short.

### Minor

- **`saturating_add` was untested at its boundary** — `wrapping_add` survived the suite.
  Fixed by `non_reference_reads_saturates_rather_than_wrapping`; the mutant now fails it.
  *(reliability)*
- **`reach` was pinned at three points where a property was available.** Both fixtures
  start at 10; a defect at `start == 0` or over a long span passed all three. Fixed with a
  proptest over the well-formed domain, bounded below the ceiling — the one documented
  exception. *(reliability)*
- **`an_inverted_region_reaches_its_own_start` claimed agreement with production in prose
  and asserted a literal.** Fixed: it now calls `production_reach`. *(reliability)*
- **The `GenomeRegion::len` defect was recorded only in a neighbouring module's test
  prose**, with no marker where a reader meets it, no owner and no removal condition — and
  in release, where overflow checks are off, the result is not a panic but a length of 0,
  a region at the ceiling reporting itself empty. Fixed: `types.rs` now carries the
  warning and names the test that pins it. The arithmetic itself is untouched — a shared
  type, out of this step's blast radius. *(smells)*
- **The crate has two names for one quantity** — `non_reference_reads` here, `alt_reads`
  in the pre-pass. A blanket rename is not free, since `min_alt_obs` is production's
  carried parameter name. Fixed by naming the synonym in the doc so a grep works from
  either side. *(naming)*
- **`reach`'s doc said "locus" where the arch says "observation"**, colliding with spec
  §1.3's reserved sense of *locus*. Fixed. *(naming)*

### Nits — applied

`production_reach` (was `productions_reach`, an unspellable possessive that reads as a
plural) and the test renamed with it; the "every partial" universal narrowed, since
`from_left`/`from_right` can clamp to a run covering the whole locus and still answer
`Partial`; the bare `(spec §3)` spelled out as `locus_generation.md §3`; the
unwrap-and-rewrap in `reach` replaced by `Position`'s own `Ord`; the "Empty and absent"
test doc trimmed to the one state it exercises; `production_reach`'s position argument
read from the fixture rather than hardcoded; and the canonicalisation the predicate
silently depends on written down.

### Nits — not applied

- The `..locus(region, observed)` functional update in the new fixture. The agent filed
  it as a nit itself: `locus` spells all six fields, so a new field is a compile error
  there.
- A destructuring rebuild of that fixture — same reason.

## 5. Out of scope observations

- **`GenomeRegion::len` is wrong at `end == u64::MAX`** and two pre-existing call sites
  (`num_obs_along_locus`, `add_generic`) would still hit it. The fix is one character
  class — `end.saturating_add(1)` — in a shared type. Worth its own commit.
- **No census test reaches `add_generic`'s multi-base path**: every generic fixture there
  is single-base. Pre-existing, untouched by A2.

## 6. What's good

- The plan's safety argument for the census move **holds and was checked, not assumed**:
  four census tests assert on the sparse allele list the moved predicate fills, and one of
  them asserts both that a differing read lands and that a matching read does not, so a
  predicate stuck in either direction fails it. Inverting the moved predicate fails four.
- The move is provably behaviour-preserving at the type level: both sides are
  `Box<[u8]>`, both spellings deref to `[u8]` and resolve to the same
  `<[u8] as PartialEq<[u8]>>::eq`, and the call site's `&Box<[u8]>` has exactly one
  coercion target.
- `reach` avoiding `GenomeRegion::len` was checked and is **CHECKED-CORRECT**: the agent
  reproduced the panic (`panicked at src/ng/types.rs:94:9: attempt to add with overflow`)
  and confirmed that mutating `reach` to production's spelling kills the ceiling test with
  that same panic.
