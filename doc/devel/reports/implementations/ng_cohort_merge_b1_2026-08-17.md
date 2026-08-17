# ng cohort merge — B1: projection onto the locus span

*Implementation report, 2026-08-17. Step B1 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §4.2 and [arch](../../ng/arch/cohort_merge.md) §4.*

## 1. Plan

Widen each member's observed sequence to the full locus span, padded with the reference
bases either side (spec §4.2, *Projection*). Nothing unifies yet — two members that widen
to the same bytes are still two projections after this step; making them one allele is
B2's.

## 2. Where the padding bases come from — a departure from the spec's sentence

Spec §4.2 says the padding bases "travel on the observation already", naming
`SampleLocusObservations::reference_bases`. **They do, but only over the observation's own
region.** A SNP at one position carries one reference base; it knows nothing of the four
its neighbour's deletion covers, and those four are exactly what the projection has to pad
with. So the step needs one thing the spec's sentence does not name: **the locus's
reference, gathered across its members**, and then each member padded from that.

The gather is sound because closing guarantees the members cover the locus with no gap: a
locus opens at some observation's first base, and every observation that joins starts at or
before the reach the locus had at the time (spec §4.1), so the members' regions form one
unbroken stretch from the locus's first base to its last. **The code asserts it rather than
resting on the argument** — see §4's uncovered-locus test.

The observations remain the only source. Nothing is fetched and no reference file is read,
which is what §4.2's sentence was really protecting.

Production widens from a reference it holds directly (`project_local_allele`,
`var_calling/per_group_merger.rs:1201-1211` — prefix, allele, suffix). The arithmetic here
is that function's; only the reference's provenance differs.

## 3. Changes made

New file [`build.rs`](../../../../src/ng/run/cohort_merge/build.rs):

- **`LocusReferenceBases`** — the reference's own bases over one whole cohort locus, with
  the region they cover. `over(&ClosedLocus)` gathers them from the members; `bases()` is
  also **the reference allele**, which is what makes B2's table hold it without a special
  case.
- **`LocusReferenceBases::placing(member) -> MemberPlacement`** — where one member sits
  inside the locus, worked out once. `MemberPlacement::project_into(sequence, &mut Vec<u8>)`
  is the widening: the reference before the member's region, the sequence the reads showed,
  the reference after it. The buffer is cleared and refilled, so one buffer serves a whole
  locus (the project's scratch-over-allocation preference; B2 allocates only when a
  projection turns out to be a new allele).
  **The handle is the review's doing and it earns its place**: taking the member and the
  sequence as two loose arguments let any member of the locus be paired with any other
  member's sequence, which pads at the wrong offset and yields a well-formed allele with
  nothing to say so.
- **`MemberPlacement::projectable_sequences()`** — the member's complete sequences, which
  is the way in that cannot reach a partial (§4).
- **`offset_within(locus_region, member_region)`** — where the member starts inside the
  locus, through `checked_sub`. Were it open-coded, a member starting earlier would wrap in
  the release profile — overflow checks off — and index the gathered reference near
  `usize::MAX`.

In [`close.rs`](../../../../src/ng/run/cohort_merge/close.rs): `span_of` is now
`pub(super)`, so the width a locus is judged on and the width its members are projected
over are one number rather than two. Its doc says so, and no longer carves out an exception
for an observation's own span — that exception was an invitation to the one mutation the
suite did not catch.

In [`mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs): `pub mod build;` and the module
doc's inventory of what has landed.

**Eight release-level assertion sites, of two kinds.** Five are against a caller's mistake —
the locus's verdict, and a member on another contig, reaching past the locus, or starting
before it — and each guards a *silent* failure: a partly gathered reference pads every allele
in the locus with `NUL` bytes, and those alleles still unify, still count, and still look
like sequences. Two are about the data an observation carries: its reference width matching
its own region, and two overlapping members agreeing on the reference. The last refuses a
partial sequence. The costs are one pass over the gathered span (bounded by the verdict
assertion to `max_cohort_locus_span` for a generic locus and the catalog's tract width for an
STR one), one comparison per base of it, and three comparisons per member.

## 4. Assumptions — one, and it is inherited rather than new

**A partial observation is not projected.** A
[`Partial`](../../../../src/ng/locus_generation/witness.rs) observation's `bases` stop
where its read's witness stopped, so padding them from the locus's reference would report
that the read showed reference bases over ground it never saw — an allele no molecule
carried. Partials are not offered:
`MemberPlacement::projectable_sequences()` yields the complete ones, the same subset
`non_reference_reads` counts over, and `project_into` panics as a backstop for a caller
that reaches into `observations` itself.

**This is the line A2 already drew, not a new decision:** `non_reference_reads` counts over
the complete observations only, for the same reason, and `PROJECT_STATUS.md` already
carries it as an owner call ("`non_reference_reads` counts complete observations only, and
the design documents do not choose"). B1 extends the same line to projection and adds
nothing to the question. What a partial should contribute waits on the censored likelihood
that does not exist yet (`spec/locus_generation.md` §3).

## 5. Tests added — 21, all in `build.rs`

The plan names two; both are the first two here. **Nine of the 21 were added after the
review**, each because a mutation of the code survived the first twelve — those rows are
marked *(review)* and the mutant each one kills is in the
[fix report](../reviews/fixes_applied_ng_cohort_merge_b1_2026-08-17.md).

| test | what it pins |
|---|---|
| `a_snp_inside_a_deletions_span_projects_onto_the_whole_span` | the plan's first case. Locus `ACGTA` at 10–14; sample 0's deletion recorded `A` over all five bases, sample 1's SNP `G`→`T` at 12 alone. The SNP projects to `ACTTA`, the deletion stays `A` |
| `an_insertion_projects_at_its_anchor_base_and_leaves_the_span_alone` | the plan's second. `GTTT` at anchor 12 projects to `ACGTTTTA` — 8 bases of allele over 5 of reference — and the locus is still 5 wide, which is why an insertion cannot push a locus past `max_cohort_locus_span` (spec §3.1) |
| `the_reference_is_gathered_across_members_that_each_cover_part_of_the_locus` | **no single member covers the locus**: 10–12 and 11–14 gather to `ACGTA`, and the left member's projection takes its two right-hand bases from the *other sample's* member |
| `one_samples_two_observations_both_reach_the_gather` *(review)* | the gather's **inner** loop. The test above spreads coverage across samples; this spreads it within one, which is spec §4.2's "two of its own observations" and the case B3 sums support over |
| `a_member_that_matched_the_reference_projects_to_the_locus_reference` *(review)* | `bases()`'s claim to be the reference allele, as a round trip, at a non-zero offset |
| `a_member_covering_the_whole_locus_projects_to_its_own_bases` | padding adds nothing at a one-position locus — and the gathered reference is asserted, so the test cannot be satisfied by code that ignores the reference entirely |
| `a_reference_containing_n_bases_is_gathered_like_any_other` *(review)* | the sentinel's one property. `N` is what every assembly gap gathers, so a sentinel spelled `b'N'` would refuse a locus its members do cover |
| `a_reused_buffer_holds_only_the_latest_projection` | the `clear()` is load-bearing: a 13-base insertion allele followed by a 1-base deletion allele, so appending *or* truncating around the old contents would show |
| `a_locus_at_the_coordinate_ceiling_is_gathered_and_padded_on_its_true_width` | both the gather and the projection measure through `span_of`, not `GenomeRegion::len()`, which answers 0 at the ceiling in release. **The member ending *at* the ceiling is now projected too** — for the SNP one base short of it the two spellings agree, so the first version pinned only the gather |
| `a_locus_its_members_do_not_cover_is_refused` | a gap leaves a byte that is not a base in every allele; the message names the position |
| `a_member_starting_before_the_locus_is_refused` | the wrapping subtraction |
| `a_member_reaching_past_the_locus_is_refused` | the same mistake from the other side, named rather than an out-of-range slice |
| `a_member_whose_reference_does_not_cover_its_region_is_refused` | `SampleLocusObservations` has public fields and nothing upstream enforces this; the consequence is a plausible allele padded from the wrong bases |
| `a_member_on_another_contig_is_refused` | every position exists on every contig |
| `two_members_disagreeing_on_the_reference_where_they_overlap_are_refused` *(review)* | samples called against different references; without the check the gathered reference is a mixture and every allele is padded from it |
| `a_failed_locus_is_not_gathered` *(review)* | spec §3.2 in code: closing is uncapped, so a locus the caller refused could allocate thousands of bases for a reference nothing will use |
| `projecting_a_member_from_another_contig_is_refused` *(review)* | the **projection's own** contig guard. Removing it was the one mutation that produced no panic at all — a member from another contig came back padded from this locus's reference |
| `projecting_a_member_reaching_past_the_locus_is_refused` *(review)* | the projection's reach guard, whose message is otherwise a raw out-of-range slice |
| `projecting_a_member_starting_before_the_locus_is_refused` *(review)* | the wrapping subtraction on the path a caller reaches directly |
| `a_partial_sequence_is_not_offered_for_projection` *(review)* | §4's line, as the API rather than as a panic: a member carrying one complete and one partial sequence offers only the complete one |
| `projecting_a_partial_sequence_is_refused` | the backstop behind that iterator |

**The fixtures spell their bases out**, unlike `close.rs`'s, whose filler `A`/`C` sequences
suit a walk that reads positions and counts and never bases. The fixture builder checks its
own reference width against its region, so a fixture a base short fails in the fixture
rather than inside the code under test.

## 6. Validation

Run in the container (`./scripts/dev.sh`), on this branch:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean
  (`Finished dev profile`). **`--all-targets` is red on this branch and was red before any
  of this work** — 49 pre-existing errors in `examples/`, `benches/` and other modules'
  test code, two of them in `src/ssr/`, which ng may not edit. The standing item under the
  block's `Open:`.
- `cargo test --lib ng::run::cohort_merge` — 59 passed, 0 failed (38 before this step; 50
  when the review started).
- `cargo test --lib` — see the commit; the whole library suite green.
- **Every mutation the review found surviving now fails exactly one test**, re-run one at a
  time from a pristine copy — the table is in the
  [fix report](../reviews/fixes_applied_ng_cohort_merge_b1_2026-08-17.md).

## 7. Tradeoffs and follow-ups

- **One allocation per locus for the gathered reference**, plus one growable buffer per
  caller for the projections. The spec's memory section (§8) prices the observation cache,
  not this; at a generic locus the gather is at most `max_cohort_locus_span` bases.
- **Overlapping members that disagree on the reference are refused** — the first draft let
  the last one copied win, and the review's point stands that this is the same
  plausible-allele failure the other checks exist to catch. It costs one comparison per base
  of the gathered span. Where the *run-level* fault belongs — samples called against
  different references — is still not this step's to decide; what changed is that it can no
  longer pass through silently.
- **Unification is B2's**, including the deletion written at two placements. B1 leaves two
  identical projections as two.
