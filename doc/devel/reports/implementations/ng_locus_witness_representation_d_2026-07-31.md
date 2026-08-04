# ng — the locus witness representation, Milestone D (consumers and surfaces)

*Implementation report, 2026-07-31. Plan:
[locus_witness_representation.md](../../ng/impl_plan/locus_witness_representation.md) Milestone D.
Design: [spec](../../ng/spec/locus_witness_representation.md) §1, §4, §8;
[arch](../../ng/arch/locus_witness_representation.md) §1.1, §2. Branch `ng-pileup-generator`,
worktree `pop_var_caller-ng-pileup`.*

**Status: Milestone D complete — D3, D4, D5, D6.** This report is extended per step and committed
with each of them. D1 and D2 landed inside C2, which could not compile without them.

**The baseline this milestone starts from**, re-measured rather than inherited: `cargo fmt
--check` and `cargo clippy --all-targets --all-features -- -D warnings` clean, `cargo test --lib
--bins --tests --examples --all-features` **2,835 passed / 0 failed**, `ng::locus_generation`
**304 passed**. The STR oracle is `tmp/witness_baseline/ssr_dump_outside_tract.tsv` (8,138 lines,
the C0 rebaseline), **not** `ssr_dump_a2.tsv`.

---

## D3 — the constructor set, reshaped by an owner decision

### What the step was asked for, and what it is

The plan asked for `ReadWitness::from_run(offset_in_locus, positions_covered, locus_len)` — the
interior-run constructor the deferred note on the variant asked for, since `from_left` and
`from_right` can only place a run *against* a border and neither can express one touching
neither.

It was raised with the owner before any code, because the step carried two questions the spec
and arch did not settle (the Milestone C structure review, F10, flagged both for exactly here):
whether arch §1.1's "all constructors return `Complete` when the clamped run covers the whole
locus" should be implemented, and whether "flush at both borders is not the same as pinned"
needed anything beyond the prose already on the type. The conversation that followed replaced
the step.

**What landed: `ReadWitness::from_witnessed_runs(runs, locus_len) -> Option<ReadWitness>`**, and
a rule that splits the constructors by **what the caller is claiming**:

- **A reach** — `from_left` / `from_right`: *the read got at least this far from this border*. A
  lower bound. On the STR path it is counted in **read** bases (`ssr.rs:897`, `tract.end -
  tract.start`) against a locus measured in **reference** positions (`ssr.rs:898`,
  `locus.segment.tract_len()`), and stutter makes those two rulers diverge, so a reach at or past
  the locus length says the read *ran out of read* — not that it reached the far border. These
  never answer `Complete`.
- **A witnessed set** — `from_witnessed_runs`: *these are the positions the read witnessed*, on
  the locus's own ruler. Completeness is then arithmetic rather than inference, and it is decided
  on the **total** positions covered, never on the outer edges.

`from_run` is **not built**: an interior run is `from_witnessed_runs([(3, 7)], len)`, and two
spellings of one run differing only in whether they decide completeness is a coin-flip for the
caller. `Complete` stays a **bare** variant that a caller writes when it knows structurally — the
STR delimiter reporting both borders of the tract anchored in this read (`ssr.rs:838`).

### Why the arch's contract was replaced rather than implemented

The implementation report for C justified the departure on *output movement*: under the contract,
`ssr::tally::tests::an_expanded_allele_merges_the_two_sides_into_one_observation` (`ssr.rs:1403`,
`from_left(9, LocusLen(6))`) turns two partials into two completes, which moves the STR dump's
`obs_complete=`/`obs_partial=` header line **and** the `depth` column on every row of the locus.

That is true but it is the weaker argument. The stronger one is correctness: `Complete` is
defined as "the read reached **both** borders and witnessed every position between them", and it
is the gate on `complete_observations()` — what a likelihood may score as an **exact** allele
length. A read anchored at one border whose read-coordinate reach happens to equal the reference
tract length has not reached the second border. Implementing the contract would score a lower
bound as a measurement.

### Why `Complete` gained no payload

The owner asked whether `Complete` should carry the locus span, so every predicate could derive
uniformly from stored data instead of branching on the variant. Three reasons it did not:

1. **The type refuses to store a locus length, deliberately.** Its own recorded reason
   (`witness.rs`, the `Partial` variant): a run clamped against *some* `LocusLen` proves nothing
   about the locus it is finally attached to — `ReadWitness` cannot know its own locus — so the
   real check lives in `num_obs_along_locus`, where the region is in hand. A stored span is
   exactly such an unverifiable claim, and it would make `Complete { locus_len: 6 }` on a
   10-position locus expressible where today `Complete` cannot be wrong: it is true *relative to
   whatever locus it is attached to*.
2. **It is paid 1.6 million times to help 872 observations.** 1,646,289 of 1,647,161 observations
   on the chr1 run are `Complete` (spec §3.1). `read_witness` is part of `ObservationKey`, the
   identity that decides which reads merge into one observation, so each would build a run and
   then be **compared and hashed as a slice** where today it is a discriminant check.
3. **Spec §3.1 settled it** — `Complete` stays a variant so `complete_observations` is a cheap
   equality — so changing it is a spec edit, not a plan note, plus the 67 sites naming the variant.

A `from_complete()` constructor was also considered and skipped: with the variant public at 69
sites, it would be a second spelling of the same value that nothing can stop drifting.

### What the code does

`from_witnessed_runs` clamps each run into `locus_len`, drops a run left covering nothing (so a
caller reaching past the locus builds a *shorter* witness rather than one claiming positions the
locus does not have, and one out-of-locus run does not sink the whole set), canonicalises through
`WitnessedLocusPositions::from_half_open_runs`, and answers `Complete` when
`positions_covered() == locus_len`. `None` when nothing survives.

`witness_of` (`open_record.rs`) keeps the part that is genuinely the fold's — intersecting each
run with the final footprint, rebasing reference positions onto the locus, the `u32 → u16`
narrowing and its panic message — and **delegates the completeness decision**, so the rule has one
home. Its four trailing code lines became one call.

**The `LocusLen` it passes is honest, and it is not the C review's F3 coming back.** F3 was about
narrowing run *offsets* through a type that means "a locus length". Here the quantity genuinely is
one: `finalise` emits the region `record_pos ..= record_end_exclusive - 1`
(`open_record.rs:836-843`), whose `len()` is `end + 1 - start` (`types.rs:93-95`) — exactly
`record_end_exclusive - record_pos`. The width of a finalised footprint *is* that locus's length.

### Departures from the plan, recorded

One, and it is the step itself: `from_run` was replaced by `from_witnessed_runs`, by owner
decision on 2026-07-31, together with the rewrite of arch §1.1's contract. Recorded in the plan's
D3 line and in the arch, the way C0 was.

### How we know it works

**Four mutations, each failing the test that names it** — run in one container start, each
applied, tested, and reverted:

| mutation | what it produced | failed |
|---|---|---|
| decide `Complete` on `span()` | the spliced read `[(0,3), (7,10)]` on a 10-position locus declared complete | `a_set_covering_every_position_is_complete_and_a_hole_is_not`, **plus 6 more: 5 existing `witness_of` fixtures and the parity anchor `every_divergence_from_production_is_one_of_the_six_named_classes`** |
| decide `Complete` on flushness | same | same 7 |
| drop the run-end clamp | `[(8,40)]` stored verbatim; and `[(12,20)]` on a 10-position locus became **`Some(Complete)`** | `runs_are_clamped_into_the_locus_and_empty_ones_dropped`, `a_set_with_nothing_inside_the_locus_answers_none` |
| drop the empty-run filter | one out-of-locus run sank the whole set to `None` | `runs_are_clamped_into_the_locus_and_empty_ones_dropped` |

The last two are caught **only** by the new tests, which is what earns them their place. The first
two are caught by the fold's fixtures as well, which is the defence in depth the delegation buys.

The holed case is the discriminating one in the first two, and it is drawn from the change's whole
purpose: a witness flush at **both** borders that covers 6 of 10 positions. The whole-locus case
alone passes under both mutations.

**The oracles.** The STR dump on tomato `SRR7279503` chr01 is **byte-identical** to
`ssr_dump_outside_tract.tsv` — 8,138 lines, zero diff — as it must be, since no STR call site
moved. `parity::ng_agrees_with_production_where_production_fabricated_nothing`,
`ng_emits_the_same_bytes_in_a_second_process` and
`every_divergence_from_production_is_one_of_the_six_named_classes` are green.

**Counts:** `ng::locus_generation` 304 → **308**; the suite 2,835 → **2,839**. `cargo fmt --check`
and `cargo clippy --all-targets --all-features -- -D warnings` clean.

---

## D4 — the surfaces: one derivation, four spellings decided, and a guard that could not fail

### What was already done, and what was left

The plan asks for two things here, and C2 had already delivered part of one. The generic dump
**already** printed one `<offset>+<positions>` per run (`partial:0+2,4+6` renders a hole as the
two runs it is, not as the span that swallows it) and **already** checked its invariant per run.
So D4's work was: the tag, the label drift, and — the part nobody had — making the per-run check
falsifiable.

### The retired variant name, finally out of the output

`observed:<offset>+<positions>` became `partial:<offset>+<positions>`, and the counters
`rows_observed` / `reads_observed` became `rows_partial` / `reads_partial` (9 sites, including the
header line the dump prints and the module doc's accounting identity). These were the last
user-visible uses of the name spec §3.1 retired: next to `complete`, "observed" is not a contrast,
because a complete witness was observed too.

### The label drift, decided rather than inherited

Three STR dumps carried the side derivation with a **byte-identical seven-line comment** and had
drifted: `ng_ssr_loci_dump` emitted `partial:left` / `partial:right`, while
`ng_ssr_cohort_stutter` and `ng_ssr_aligner_bakeoff` emitted `partial_left` / `partial_right` —
and all three said `partial:interior`, so two of them mixed both separators **inside one
function**. A consumer grepping `partial:` got one tool's sides and another tool's interiors.

The derivation now lives once, in `examples/shared/witness_side.rs`, reached by
`#[path = "shared/witness_side.rs"] mod witness_side;` in each of the three (cargo discovers
`examples/*.rs` and `examples/*/main.rs`, so a plain file in a subdirectory is compiled only where
an example asks for it). It returns a `WitnessSide` enum; **each tool keeps its own strings**,
because a dump's output is its own contract and must not move because a sibling's did — which is
the plan's instruction, and it is what let the drift be *decided*: every tool now spells the colon
form, the one that was already internally consistent.

**Exactly which bytes moved.** `ng_ssr_loci_dump` — **none**; it already spelled the colon form,
which is why the byte-identity oracle is intact. `ng_ssr_cohort_stutter` and
`ng_ssr_aligner_bakeoff` — two string literals each, in the `coverage` column: `partial_left` →
`partial:left`, `partial_right` → `partial:right`. `complete`, `partial:interior`, `no_border` and
`capped` unchanged; no other column, count or header touched. Neither tool has tests or a
committed baseline, so there is nothing to rebaseline — but they do have a downstream consumer.

**The consumer, found by grepping for the old strings rather than assumed absent.**
`benchmarks/ssr_hg002/scripts/ng_ssr_aligner_bakeoff_dashboard.py` maps the `coverage` column into
outcome classes, and its map was keyed on `partial_left` / `partial_right`. An unmapped label
becomes `NaN`, which every downstream count silently drops — so the rename would have cost the
notebook its whole *partial* class with no error and no visibly wrong plot. The map now carries
both spellings (old TSVs still classify) plus `partial:interior`, and an assertion turns an
unknown label into a named failure. The other dashboard that reads this column,
`benchmarks/ssr_tomato1/scripts/ng_ssr_cohort_stutter_dashboard.py`, filters only
`coverage == "complete"` (`:172`) and is unaffected.

### The per-run bound had nothing aimed at it

C2 made the generic dump's invariant a check per run, correctly, and **no test could tell it from
either enclosing formula** — every locus the walk produces satisfies all of them. Two tests now
build a locus by hand (`SampleLocusObservations`' fields are public, and
`WitnessedLocusPositions::from_half_open_runs` canonicalises without clamping, which is the
documented split — the clamp belongs to `num_obs_along_locus`, where the region is in hand):

- a witness of `[(0,2), (4,12)]` on a 10-position locus **must panic**. It is the case that
  discriminates: it covers 2 + 8 = 10 positions starting at 0, so the pre-C2 formula the plan
  names — `offset_in_locus + positions_covered <= footprint`, i.e. `0 + 10 <= 10` — *passes* it,
  while its second run genuinely runs two positions past the locus;
- a witness of `[(0,2), (4,10)]` passes, and renders `partial:0+2,4+6`, so the first test is
  failing on the overrun rather than on the fixture's shape.

### How we know it works

**Three mutations, each failing the test that names it:**

| mutation | result |
|---|---|
| replace the per-run loop with the enclosing formula (`first_run().0 + positions_covered()`) | `a_run_reaching_past_the_footprint_is_caught_even_when_the_totals_fit` — *"test did not panic as expected"* |
| swap `Left` and `Right` in the **shared** derivation | `ng_ssr_loci_dump`'s `the_fixtures_partials_are_asymmetric_and_so_can_catch_a_side_swap` fails (`left: 1, right: 2`) — so the one body is under test on behalf of all three tools, which is the second thing sharing it bought |
| put `observed:` back on the generic dump | 6 tests fail (of 13) |

**The oracles.** The STR dump on tomato `SRR7279503` chr01 is **byte-identical** to
`ssr_dump_outside_tract.tsv` (8,138 lines, zero diff) — the label unification landed on that
tool's own spelling, so the oracle did not move. The three parity anchors are green.

**Counts:** the suite 2,839 → **2,841**; `ng::locus_generation` unchanged at 308 (D4 is all
example-side). `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings`
clean.

---

## D5 — the census counts the holes, and is honest about reading zero

### What the step is

`DivergenceCensus` gains `holed_witness_reads` and `hole_positions`: the reads whose witness is
more than one run, and the locus positions inside those runs' gaps. Both are weighted by
`num_obs`, like every other census deliverable, because an observation is shared by the reads that
agreed.

The positions come out as **`span() - positions_covered()`** — the distance the runs sit in, minus
what they cover. That pair of accessors exists precisely so this is not open-coded as "last end
minus first start", which reads as a coverage and is not one (Milestone C review, F7). The
difference is non-zero exactly when the witness has more than one run: the set is canonical, so
two runs are separated by at least one unwitnessed position, and one run has nothing between
anything.

Both are printed by both census reports — the synthetic differential's and the real-data one.

### The one design decision: no floor, a positive control instead

Every other census deliverable is floored — `fabricated_reads > 0`, `stale_widen_reads > 0` — on
the argument that a measurement read off production's observations rather than off the
classification can silently stop measuring. **These two cannot be**, because their expected value
on every alignment this repo can currently run is **zero**: spec §8 measured 0 holed witnesses in
225 million DNA-seq event-folds, and structurally so — a `Skip` emits no event, so an intron
cannot widen a record on its own, and modern Illumina puts `N`s at read ends where they cannot
make a hole. The number that matters is RNA-seq's, and no spliced alignment was available (spec
§8, still open; plan E4).

That is exactly the shape spec §8 warned about in as many words — *"zero is also what a miswired
probe reports"* — so the guard is a fixture that produces the thing the counters count:

- a read blind over three positions in the middle, witnessing `[(0,3), (6,10)]` of a 10-position
  locus, carried by 4 reads → **4 holed reads, 12 hole positions**;
- a one-run partial beside it → **neither counter moves**, while `fabricated_reads` and
  `fabricated_ref_bases` still do (7 reads, 42 bases). Without that half, "holed" could collapse
  into a second spelling of "partial" and the counter would be measuring a class that already had
  a name.

### How we know it works

| mutation | what it produced | test result |
|---|---|---|
| measure the span rather than the gap (`span() - span()`) | the counters report nothing at all | `the_census_counts_a_hole_and_the_positions_inside_it` — `left: (0, 0)`, `right: (4, 12)` |
| count every partial as holed (`if true`) | the one-run partial counted | `left: (7, 0)`, `right: (0, 0)` |
| drop the per-read weighting (`+= 1`) | observations counted instead of reads | `left: (1, 12)`, `right: (4, 12)` |

**The STR oracle was not re-run for this step, and could not have moved:** `mod parity` is
`#[cfg(test)]` (`pileup/mod.rs:131-133`), so none of it is compiled into the release example the
oracle runs. It is re-run once at the end of the milestone.

**Counts:** `ng::locus_generation` 308 → **309**; the suite 2,841 → **2,842**. `cargo fmt --check`
and `cargo clippy --all-targets --all-features -- -D warnings` clean.

---

## D6 — the spliced fixture, permanent

### What the step is

The regression anchor for the whole change, and the only fixture in this milestone drawn from a
real failure rather than constructed to exercise a branch. Two reads over a 60-base reference:

- **the spliced read**, `3M 15N 3M` from 28 — exon 1 at 28–30, a 15-base intron at 31–45, exon 2
  at 46–48. A `Skip` emits no event, so the read witnesses six positions in two runs;
- **the deletion read**, `3M 20D 3M` from 26 — it anchors a record at 28, the base before the
  deletion, and widens its footprint to `28..=48`, twenty-one positions.

The deletion read is why the hole is *inside a record at all*: an intron cannot widen a record on
its own, so without an indel allele spanning it the two exons would simply be separate records and
there would be nothing to represent (spec §8).

**The assertions** are the two halves of the failure. The spliced read's observation is *there* —
before C3 it was absent from the record entirely, `apply_events_into` having answered `None` for a
non-contiguous witness — and its witness reads `[(0,3), (18,21)]`, two runs, rather than a span
that swallows the fifteen positions it never saw. Its bases are the two exons and nothing from the
intron. `reads_without_observation` is 0: a holed witness is evidence, not a read that witnessed
nothing.

### The knife-edge, asserted from both sides

The plan asks for the one-position sensitivity to be *recorded in a comment*. It is asserted
instead, because a comment about a boundary is not a test of one:

- a **17**-base deletion widens the record to `28..=45`, one position short of exon 2, and the
  spliced read's witness inside it is exon 1 alone — one run, which the old representation
  described perfectly;
- **18** reaches 46, the first base of exon 2, and the same read is holed at once.

So the change earns its keep on a single deleted base. *(The plan says the boundary is at 16; that
was the throwaway probe's geometry. This fixture's own numbers are 17 and 18, verified by running
it, and they are what the test and the plan now record.)*

### How we know it works

**The mutation is C3 put back** — the fold reporting "nothing witnessed" whenever the witness is
more than one run, which is exactly the pre-C3 behaviour:

```
test …pileup::tests::a_spliced_read_across_a_widened_record_is_recorded_with_both_of_its_runs ... FAILED
test …pileup::tests::one_more_deleted_base_is_what_turns_the_spliced_read_into_a_holed_one ... FAILED
test …open_record::tests::a_read_with_a_hole_is_counted_neither_as_capped_nor_as_witnessing_nothing ... FAILED
test …open_record::tests::a_read_whose_witness_splits_when_the_record_widens_stays_in_it ... FAILED
test …open_record::tests::a_read_folding_at_four_positions_of_one_record_is_one_observation ... FAILED
test result: FAILED. 304 passed; 7 failed
```

Both new fixtures fail, alongside C3's five. And D3's mutations already showed that deciding
completeness on the outer edges or on `span()` turns this same shape into a `Complete` — so the
fixture is pinned against both ways of losing a hole: dropping the read, and describing it as
whole.

**The oracles, at the end of the milestone.** The STR dump on tomato `SRR7279503` chr01 is
**byte-identical** to `ssr_dump_outside_tract.tsv` — 8,138 lines, zero diff — across all four
steps. `parity::ng_agrees_with_production_where_production_fabricated_nothing`,
`ng_emits_the_same_bytes_in_a_second_process` and
`every_divergence_from_production_is_one_of_the_six_named_classes` are green.

**Counts:** `ng::locus_generation` 309 → **311**; the suite 2,842 → **2,844**. `cargo fmt --check`
and `cargo clippy --all-targets --all-features -- -D warnings` clean.

---

## D7 — the hole counters where a real BAM can reach them (added at Checkpoint D)

### What the step is, and why it was not already done

D5 put the two hole counts on the divergence census. The Milestone D structure review pointed out
what that costs and the owner called it: the census is `#[cfg(test)]`, and it only measures loci
where **production's** walker also produced a record. So it can never answer the one question this
whole representation was built to answer — how often a read sees a locus in two pieces on real
RNA-seq (spec §8, open; plan E4).

The same two counts now ride the walk's own `RunSummary` → `PileupGeneratorCounts`, and
`ng_generic_loci_dump` prints them in its header line. Pointing that tool at a spliced BAM answers
the question, with no probe and no comparison against the old caller.

### Where they are counted

In `finalise`, in the loop that already resolves each folded read's witness against the final
footprint — so the walk pays one subtraction per partial read and nothing per complete one. A read
is holed when `span() - positions_covered()` is non-zero, which is exactly "more than one run",
since canonical runs are separated by at least one unwitnessed position.

Two exhaustive destructures made the wiring self-checking rather than hopeful: adding the fields
to `RunSummary` failed the build at `PileupGeneratorCounts::fold_region_walk` and at `parity.rs`'s
`ng_counters`, which is where the decision "production has no counterpart, so bind and drop by
name" had to be made explicitly. That is the third counter to take that route.

### How we know it works

The end-to-end test is **D6's own fixture through the real walk**, and its numbers are exact
rather than "greater than zero":

- at a 20-base deletion — the record widened across the intron — **1 holed read, 15 blind
  positions**, the intron being 31–45;
- at a 17-base deletion — the footprint stopping one position before exon 2 — **0 and 0**.

The second half is what makes the first discriminating: a counter that simply counted every
partial read would pass the first assertion and fail the second.

**Counts:** `ng::locus_generation` 312 → **313**; the suite 2,847 → **2,848**. STR dump
byte-identical to the C0 baseline; `cargo fmt --check` and
`cargo clippy --all-targets --all-features -- -D warnings` clean.

### And a measurement the owner's question asked for

`ReadWitness` stores locus positions — reference coordinates relative to the locus start — so
nothing is held in read coordinates. Read coordinates enter only as the repeat's *length in the
read*, which the delimiter returns (`RepeatSpan`, a read-space range) and which the STR path turns
into how far into the tract the read reached.

On chr01 of tomato `SRR7279503`, **2,530 of 6,216 partial rows (41%)** are reads whose repeat, in
read bases, reached or passed the reference tract's length. Their stored witness covers the tract
end to end and is honest — the read did see every reference position — and they are correctly not
`Complete`, because the read ran out and so did not measure the allele. What collapses is the
printed label: all 2,530 render as a left-edge partial, including the reads anchored on the right.
Recorded in spec §8; the decision on a fourth label is the owner's.

---

## D8 — a fourth label, and the convention that justifies the clamp (added at Checkpoint D)

### The question that produced it

The owner asked why a witness is not simply stored in reference coordinates, given the caller
always has an alignment. The answer is that it already is — a witness is a set of locus positions,
which are reference positions counted from the locus start. Read coordinates enter at one place
only, as a *length*: the tract delimiter returns the repeat's span in the read's own sequence
(`RepeatSpan`, a read-space range), and the repeat caller turns that length into how far into the
tract the read reached.

The owner's follow-up settled the part that had been recorded as unresolvable. I had written that
placing a read's extra repeat copies on the reference is arbitrary; it is not, once the anchored
side is allowed to choose the direction.

### The convention, now written down

**Lay the read's repeat down from the border it anchored.** A read that held the left flank starts
its repeat at the tract's left border; one that held the right flank ends its repeat at the right
border. Bases past the far border are extra copies — an insertion — rather than positions of the
tract. It is the same rule indel left-alignment uses, with the anchored side choosing the
direction.

Two consequences, and **the existing clamp is exactly this rule implemented**:

1. a reach shorter than the tract covers that many positions from the anchored border;
2. a reach at or past the tract length covers the tract **end to end**, with the surplus outside as
   inserted copies.

So this step changes no behaviour. What it changes is the justification: the clamp reads as a
placement rule rather than as a saturation artefact. Written where the repeat caller computes the
reach, and on `from_left` / `from_right`.

### The label, and why case 2 needed one

Covering every position is still not a measurement. The read anchored **one** border and then ran
out, so the allele can continue past what it showed and the evidence stays a lower bound. The
representation had this right — those reads are `Partial`, never `Complete`. What did not have it
right was the label, which asked "does it touch the left border?" first and answered *left* for
every one of them, **including the reads anchored on the right**.

`WitnessSide` gained `BothBorders` and the border match lost its `(true, _)` wildcard — that
wildcard is what had been swallowing the case. All three STR dumps spell it `partial:both`.

**It is not a corner:** on chromosome 1 of tomato `SRR7279503`, **2,530 of 6,216 partial
observations, 41 %**, are reads whose repeat reached or passed the reference tract's length.

### The rebaseline, verified rather than asserted

The STR dump moved on purpose, for the second time in this plan. Compared line by line against the
previous baseline rather than by `diff` block counts:

```
lines differing at the same line number: 2530
which column indexes differ: {(6,): 2530}
label transitions: {('partial:left', 'partial:both'): 2530}
```

Every changed row differs in the `read_witness` column **and nothing else**; the two header lines,
the `depth` column, the row order and every other field are unchanged. New baseline
`tmp/witness_baseline/ssr_dump_partial_both.tsv`.

*(A plain `diff | grep -c '^<'` reports 2,536, which is block alignment around the changed rows and
not six extra changes — the line-by-line comparison above is the one to trust.)*

**The downstream consumer follows:** the bake-off dashboard maps `partial:both` into its partial
class. It would otherwise have failed loudly rather than silently, on the unknown-label assertion
added during the Milestone D review.

**Counts:** the label tests in all three STR dumps gained the new case; the suite is **2,848**,
`ng::locus_generation` **313**, both unchanged by this step since the new assertions live in
existing tests.

---
---

# Milestone E — verification at scale

*Same report, continued. The milestone asks four questions of the finished change: did the STR
path move only where we decided it should, does the generic path still agree with the old caller
on real data, does the representation cost an allocation per observation, and how often does the
hole actually fire.*

## E1 — the STR path across the whole change: every difference accounted for

**The claim to test.** Spec §7.1 asks for byte-identity, with two deliberate exceptions recorded
along the way. So the real question is not "is it identical" — it is not — but **"is every
difference one of the two decisions, and nothing else?"**

**How it was checked.** The dump from the start of the plan and the dump from its end, on
chromosome 1 of tomato `SRR7279503`, aligned row by row with `difflib.SequenceMatcher` over the
rows with the label column blanked — so a relabelled row aligns to itself instead of looking like
a delete plus an insert. Then the deletions were classified, and the labels compared within the
aligned pairs.

```
observation rows: 11315 before, 8135 after  (delta -3180)

ignoring the label column: 8135 rows identical, 3180 deleted, 0 inserted, 0 replace-blocks
  of the 3180 deleted rows: 3180 carry NO observed bases; labels {'partial:left': 1608, 'partial:right': 1572}

label transitions among surviving rows: {('partial:left', 'partial:both'): 2530}  (total 2530)
```

**The account, in full.** Three header/column lines and two classes of row change, and that is
everything:

| what moved | how much | the decision |
|---|---|---|
| the column header `read_coverage` → `read_witness` | 1 line | the rename (Milestone A) |
| `reads_without_observation` 2,561 → 9,265, `obs_partial` 13,789 → 7,085 | 2 header fields | C0 — reads that clip the window and never enter the tract are not STR partials |
| rows deleted | **3,180**, and **every one carries no observed bases** | C0, the same decision |
| rows relabelled | **2,530**, all `partial:left` → `partial:both` | D8 — a witness touching both borders is not a prefix |

**Zero rows appeared, zero rows were reordered, and no surviving row differs in any column other
than the label.** `obs_complete` is 15,404 before and after; the locus count and the
zero-coverage count are unchanged.

That is a stronger statement than the plan asked for. "Byte-identical apart from X" leaves the
reader to trust that X is all there is; this enumerates every byte that moved and names the
decision behind it. The two moves are also *disjoint in kind* — one deletes rows with empty
bases, the other rewrites one column of rows that survive — so neither could have hidden the
other.

**The committed fixture**, the other half of E1, is the ten tests in `ng_ssr_loci_dump`, which
run in the suite on every gate and assert the tool's output against a fixture in the tree.

## E2 — the generic path against the old caller, on real reads

The differential runs both walkers over one prepared read stream and compares record for record.
Three runs, two organisms, two depths:

| run | loci | reads | class 1 (partial witness) | the deliverable | holed reads |
|---|---|---|---|---|---|
| HG002 30×, chr1:1–6 Mb | 47,752 | 5,683 | 12 | 16 reads / 12 loci (0.03 %) / 22 bases | **0** |
| HG002 300×, chr1:1–6 Mb | 48,905 | 55,054 | 162 | **871 reads / 162 loci (0.33 %) / 1,550 bases** | **0** |
| tomato SRR7279503, chr01:1–6 Mb | 96,253 | 6,146 | 29 | 53 reads / 29 loci (0.03 %) / 393 bases | **0** |

All three green: *"every region, reference sequence and counter identical"*, and every divergence
falls in one of the six named classes.

**The 300× line is the result.** `871 reads over 162 loci (0.33 %) with 1,550 reference bases` is,
digit for digit, what the generic generator's own Milestone D reported **before any of this change
existed**. The whole witness representation — the renames, the two set types, the fold speaking in
sets, a holed read being recorded instead of discarded, the constructor split, the relabelling —
reproduces the prior measurement exactly on real high-depth data. That is the strongest available
statement that C3 changed nothing on DNA-seq, and it is a comparison against a number nobody could
have tuned toward, because it was published first.

**The hole class is counted and reads zero on all three**, which is what spec §8 predicts
structurally: a ref-skip emits no event, so an intron cannot widen a record on its own, and modern
Illumina puts `N`s at read ends where they cannot make a hole. Zero here is the *prediction*
confirmed, not a missing measurement — the same counter reads 400 reads / 528 positions on the
synthetic corpus, which does emit ref-skips, and that is what the floor asserts on.

## E3 — the memory question, and a finding

### The requirement is met

Spec §5's requirement is **"no allocation on an observation that witnessed one run"**. That is a
property of the encoding, and it is now pinned where it is decided:
`a_witness_of_one_or_two_runs_holds_them_inline_and_allocates_nothing` asserts that a one-run and
a two-run set do not spill to the heap, and that a three-run set does — so the boundary is stated,
and shrinking the inline capacity to one becomes a test failure rather than a silent cost per
observation.

### The chromosome-scale numbers reproduce

`ng_generic_loci_dump` over chr1 of HG002 30×:

```
# generic_loci=1541788 …
# rows_complete=1646289 rows_partial=872 …
# reads_with_holed_witness=0 hole_positions=0
```

**1,541,788 loci and 1,646,289 + 872 = 1,647,161 observations** — both exactly Milestone D's
numbers, and 1,646,289 complete is exactly spec §3.1's figure. Wall time is unchanged (41 s in
both trees, same machine, same container).

### ⚠ The finding: peak resident memory grew by about 70 %, and it appears at C2

Measured like for like — same machine, same container, same input, same method (`VmHWM` from
`/proc`, polled), output to `/dev/null`, and each point run twice:

| tree | what it is | peak RSS |
|---|---|---|
| `11de107` (B3) | the set types exist and nothing uses them | **501 MB**, repeat **523 MB** |
| `761d53e` (C2) | `ReadWitness::Partial` carries a set | **892 MB**, repeat **876 MB** |
| `93c7461` (D8, head) | end of the plan | **859 MB** |

**+350 MB, about +68 %, and the whole of it appears at C2** — the commit that replaced the
variant's two `u16`s with a witnessed set. The RSS trajectory is linear and monotonic through the
walk in both trees, with **no late spike**, and the slope roughly doubles: ~11 MB/s at B3 against
~20 MB/s at C2, over the same 1,647,161 rows in the same wall time. So the cost is **per row**,
and it is about 240 bytes per observation.

**What it is not.** It is not the dump tool changing: the tool's diff across C2 is the match arm
and the per-run assertions, with no new per-row allocation. It is not more rows: the row counts
are identical to the digit. It is not the row struct's shape: its fields are unchanged across C2
and its size is the same 152 bytes.

**What I have not established is the mechanism**, and I am not going to guess at it. The honest
statement is: a reproducible +68 % peak-RSS cost at chromosome scale, localised to one commit and
to per-row cost, with the cause unidentified.

**Why the instrument cannot answer it.** This measures `ng_generic_loci_dump`, whose peak is
dominated by its own whole-run row buffer — the Milestone D review established that and withdrew
the earlier RSS conclusion for exactly this reason. The generator's own live memory is not
measured here at all, and ng still has **no committed heap profile** (the perf review's M1). A
`dhat` harness over the generator would give live bytes and an allocation count attributed to a
call site, which is what this question needs and what `VmHWM` cannot give.

**Recommendation, for the owner at Checkpoint E:** do not accept Milestone E as complete on this
point. Either build the `dhat` harness and attribute the 240 bytes per observation, or record the
cost as accepted with its size stated. What should not happen is the number going into the
record as "measured, fine" — spec §5 asks for allocations per observation to be *rejected at
review*, and this is the review.

## E4 — the RNA-seq rate: still not run, and now a one-command answer

No spliced alignment is available in the tree or the benchmark bundles, so the open question — how
often a read witnesses a locus in two pieces on real RNA-seq — is still open, exactly as spec §8
left it.

What changed is that answering it no longer needs a probe. D7 put `reads_with_holed_witness` and
`hole_positions` on the walk's own run summary and into `ng_generic_loci_dump`'s header, so the
answer is one command over any spliced BAM. On all three DNA-seq runs above the counter reads 0,
which is both the structural prediction and the positive control that the plumbing reports what it
is given — the same counter reads 400 / 528 on the synthetic corpus.
