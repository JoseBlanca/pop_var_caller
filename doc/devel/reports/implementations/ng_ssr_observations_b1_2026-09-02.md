# ng STR observations — B1: the routing policy comes from the caller's flags

*2026-09-02. Step B1 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §2.1, §2.3 and §9](../../ng/spec/run_ssr_observations.md). Branch
`ng-ssr-observations`.*

## Plan

`call-from-alignments` asked the repeat catalog with `StrRepeatCriteria::default()`, which
*is* the floors the file is stored at. The file is built deliberately below every calling
floor so that a caller can put its own line anywhere inside that gap by filtering rather than
re-scanning — and the run did not honour it, so every row the file held became an STR locus
of the run, routed to a generator that does not exist yet.

Five flags, named as `type-regions` names them, defaulting to ng's measured calling floors:
`--min-copies`, `--min-period`, `--max-period`, `--max-str-len`, `--min-purity`. They build a
`StrRepeatCriteria` that `segments_over` passes to `catalog.genome_segments`.

## Assumptions

**The flank floor is not a flag**, per spec §2.1, and this is the one axis where that is a
requirement rather than a choice: the rows below the file's 15 bp were never written, so no
run can ask for less. It comes from `StrRepeatCriteria::from(&TypedRegionConfig)`, which is
where that reasoning already lived.

**The score floor is not a flag either**, and the spec does not list it. It gates the
*scanner*'s output and a catalog reader has no scanner; at its default of 0 it rejects nothing
Ruzzo–Tompa can emit.

**A `TypedRegionConfig` is built only to be converted.** Its scan half is unread — the
conversion takes the classification rules and the satellite cap and nothing else. The
alternative was to spell `min_flank_bp: Bp(CATALOG_MIN_FLANK_BP)` here, which is the same
policy in a second place; a comment says the scan half is dead rather than leaving a reader
to wonder.

**The five flags get their own help heading, "What counts as a repeat"**, where the plan says
only that they are spelled as `type-regions` spells them (`type-regions` files them under
"Advanced"). Five flags that between them define what an STR locus *is* are worth grouping.

## Changes made

| file | change |
|---|---|
| `call_from_alignments.rs` | the five flags on `CallFromAlignmentsArgs`; `routing_criteria`; `catalog_error_naming_the_flag`; `segments_over` uses both; two error variants |
| `call_from_alignments/tests.rs` | five tests, and a reference fixture with a tract on each side of the calling floor |

**Two error variants, and each exists because the alternative is unhelpful rather than
wrong.** `PeriodRange` catches a range the wrong way round — clap bounds each end to 1..=6 on
its own, so that is the only way left to type one. `RoutingBelowCatalog` names the flag that
made an unservable request: the catalog's own refusal reads *"period 1: catalog holds tracts
of 5 copies and up, reader asked for 3"*, which is two numbers and no knob. The mapping is
exhaustive over `CriteriaRefusal`, so a new bounded axis is a compile error here rather than
silently inheriting the no-flag answer. The flank axis maps to no flag and falls through to
the general catalog error, which already says the file has to be rebuilt.

## Tests added

- `the_default_routing_is_the_calling_floors_and_not_the_catalogs` — parsed through clap, not
  read off the struct, so a `default_value` drifting from the library constant fails. Asserts
  the gap itself: at every period the run now needs strictly more copies than the file holds
  from, and the satellite cap is 100 bp against the file's 500.
- `every_routing_flag_reaches_the_criteria_the_catalog_is_asked_with` — all five moved at once
  to values nothing else produces, each read back on its own. A flag that parsed and went
  nowhere is the silent failure this step was isolated for.
- `a_period_range_the_wrong_way_round_is_refused_before_the_catalog_is_opened`.
- `a_request_the_catalog_cannot_serve_names_the_flag_that_made_it` — each refusal to its flag,
  including both ends of a period range and the flank axis that has none.
- `a_tract_below_the_calling_floor_becomes_generic_ground_and_one_above_it_stays_a_tract` —
  the step's own end-to-end oracle, below.

### The end-to-end fixture, and why it is shaped that way

A one-contig reference of 136 bases: `CGTGCTG` repeated as filler, a **6**-base run of `A` at
41–46, and a **10**-base run at 87–96. The filler's period is 7 — outside the 1..=6 the
scanner looks for — and it contains no `A`, so neither homopolymer can grow into it. Forty
bases of filler each side is past both the 15 bp of flank the file requires and the 15 bp
within which two tracts would be bundled instead of being loci. The catalog is built beside it
exactly as `repeat-catalog` would build one.

Six and ten straddle the gap: the file holds period-1 tracts from 5 copies, ng calls from 8.
One switch, two settings, the same reference and the same catalog file — at the file's floors
both runs are repeat tracts; at ng's, the six-base run is `Generic` and the ten-base run is
still a tract. The fixture's own half is asserted first (both are tracts at the file's
floors), because without it the test would pass over a reference with no repeat in it at all.

**Mutation-tested**: replacing `routing_criteria`'s body with `StrRepeatCriteria::default()`
— exactly the behaviour this step removes — fails three of the five tests, the end-to-end one
among them. The mutation was reverted from a backup and its absence checked by grep before
anything was committed.

## Validation

In the dev container:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --tests --examples --all-features --no-fail-fast` — **5,920 passed, 0
  failed, 14 ignored** in the library suite (5,915 before this step); every integration target
  green.
- Still red, and unchanged by this step: three tests in the two locus dumps, recorded in
  `PROJECT_STATUS.md` — a fixture asking for a 30 bp flank against a 15 bp bundle threshold, a
  chain id on a whole-footprint reference match, and rows whose chain ids depend on where the
  analysed regions were cut. And `benches/psp_writer_perf.rs`, which indexes one past its
  fixture.

## Tradeoffs and follow-ups

- **Nothing here says where the frontier belongs.** The defaults are ng's measured stutter
  onsets; whether period-1 tracts should be on the repeat path at all is spec §9's open
  question, and the period × length measurement is spec §8's deferral. This step makes the
  line a knob and sets it; it does not measure it.
- **A supplied parameters file's criteria are not compared yet** — B2.
- **The recovery is not measured yet** — B4 re-runs the GIAB benchmark against the loss
  report's predicted ≈0.97 SNP / ≈0.94 indel recall.
