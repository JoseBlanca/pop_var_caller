# ng — the locus witness representation, Milestone D (consumers and surfaces)

*Implementation report, 2026-07-31. Plan:
[locus_witness_representation.md](../../ng/impl_plan/locus_witness_representation.md) Milestone D.
Design: [spec](../../ng/spec/locus_witness_representation.md) §1, §4, §8;
[arch](../../ng/arch/locus_witness_representation.md) §1.1, §2. Branch `ng-pileup-generator`,
worktree `pop_var_caller-ng-pileup`.*

**Status: D3 complete. D4–D6 not started.** This report is extended per step and committed with
each of them. D1 and D2 landed inside C2, which could not compile without them.

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
   equality — so changing it is a spec edit, not a plan note, plus 69 sites naming the variant.

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
home instead of two. Its six trailing lines became one call.

**The `LocusLen` it passes is honest, and it is not the C review's F3 coming back.** F3 was about
narrowing run *offsets* through a type that means "a locus length". Here the quantity genuinely is
one: `finalise` emits the region `record_pos ..= record_end_exclusive - 1`
(`open_record.rs:832-836`), whose `len()` is `end + 1 - start` (`types.rs:93-95`) — exactly
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
| decide `Complete` on `span()` | the spliced read `[(0,3), (7,10)]` on a 10-position locus declared complete | `a_set_covering_every_position_is_complete_and_a_hole_is_not`, **plus 5 existing `witness_of` fixtures** |
| decide `Complete` on flushness | same | same 6 |
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
