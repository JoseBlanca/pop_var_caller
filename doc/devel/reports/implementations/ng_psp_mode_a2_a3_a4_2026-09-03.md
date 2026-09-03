# A2+A3+A4 — the header's remaining §6.1 fields: the read-group table, the reach ceiling, the read filters

**Date:** 2026-09-03. **Plan:** [`run_driver_psp_mode.md`](../../ng/impl_plan/run_driver_psp_mode.md),
Milestone A steps A2, A3, A4 — **bundled into one loop iteration deliberately**: all three
extend the same file (`src/ng/psp/header.rs`), share the same fixtures and test surfaces, and
the milestone checkpoint follows immediately, so three near-empty loops would have re-reviewed
the same seams three times. The bundle is named here and in the commit, per the
plan-driven skill's granularity rule. **Branch:** `ng-psp-mode`, on top of A1 (`114efe24`).

---

## What landed

**A2 — the read-group table.** `Header` gains `read_groups: Vec<ReadGroupIdentity>`; each row
is the `@RG ID` verbatim, the library, and the walk-local identifier
(`walk_local_id: ReadGroupId`, renamed by the review to match the wire key) — the three
things spec §6.1 asks for and no more (the sample is
the header's own field; file paths are provenance basenames). Wire: one `[[read-group]]` row
per group with `id`, `library`, `walk-local-id`. Two decisions worth recording:

- **The identifier is stored explicitly *and* required to equal the row's position** —
  redundancy checked on both sides, like the declared length against the sentinel — so the
  number a person reads beside an id and the number code derives from order cannot disagree.
- **A duplicated `@RG ID` is legal here, deliberately.** SAM makes the id unique only within
  one file, so a sample sequenced across files may carry two rows with one id and different
  libraries. The refusal of a table calling cannot merge is the calling stage's (spec §6.2,
  plan step E1) — refusing at write would refuse legitimately-headered inputs the walk can
  process fine.
- **Rules:** the table is non-empty; identifiers are the walk's numbering from zero in order;
  ids and libraries are non-empty and hold no control characters. **Control characters, not
  all whitespace**: SAM allows a space in an `@RG` value and real archives carry them, so the
  whitespace rule the contig names use would refuse real files — pinned by
  `a_read_group_id_with_a_space_round_trips`.

**A3 — the observation reach ceiling.** `Header.observation_reach_ceiling_bp: Bp`, wire key
`observation-reach-ceiling-bp`, a top-level scalar. It is the locus generator's own span cap
(`PileupGeneratorConfig::max_record_span`), known before the first record because it is a
setting; a sizing fact, not a correctness fact — nothing consumes it until plan step E4. Rules:
non-zero (an observation covers at least one base) and within TOML's signed integer. The
minor-version test's "unknown later key" had been exactly this key; it now uses an invented
one.

**A4 — the read filters into provenance.** The enumeration lives on the config itself —
`ReadFilterConfig::provenance_parameters()` in `ng::read` (moved there by the review: the
first draft put it on `WriterProvenance`, which made the psp store import a pipeline stage) —
and psp keeps a generic seam, `WriterProvenance::record_parameters(entries)`. One key per
configurable filter: `read-filter-min-mapq`, `read-filter-min-read-length-bp`,
`read-filter-drop-qc-fail`, `read-filter-drop-duplicates`,
`read-filter-max-read-mismatch-fraction`, `read-filter-mismatch-bq-floor`; the family is
published as `READ_FILTER_PROVENANCE_KEYS` so a re-recorder can clear stale keys first.
Recorded, never compared (spec §6.1 — the census digests them). An **off** filter is the
string `"off"`, never omitted: an absent key would read as "unrecorded", a different fact.
The one deliberate omission is the mismatch base-quality floor beside an off mismatch filter,
which would invite misreading; the unconditional drops (secondary, supplementary, unmapped,
undecodable CIGAR) are not settings and so have no key. The config is destructured
exhaustively, so a filter added later fails to compile rather than going unrecorded.
`a_written_header` (the header module's own fixture) routes through the real path with
`ReadFilterConfig::default()`; the other fixtures record no filters, which is itself a legal
provenance state.

`format_version` stays (1,0) for all three — still nothing written outside tests predates them.

## Tests

- Round trips inherit the read-group table and the reach ceiling through every fixture
  (`a_written_header`, the writer's shared `a_header`, `mod.rs`, `record.rs`, both examples,
  the bench); the recorded filters ride only in `a_written_header`, through the real path.
  The parity example's exhaustive header destructure gained both fields with assertions.
- Nine rows in the two-sided rule table: zero and over-wide reach ceiling; empty read-group
  table; identifier out of walk order **and** repeating zero (the boundary case — two
  zero-numbered walks pasted together); empty `@RG ID`; newline **and** tab in an id;
  newline in a library.
- `a_read_group_id_with_a_space_round_trips` pins the deliberate narrowness of the
  control-character rule; `two_read_groups_sharing_an_rg_id_round_trip` pins the deliberate
  acceptance of duplicated `@RG ID`s (the un-mergeable refusal is the calling stage's).
- `provenance_parameters_pin_every_filter_value` (in `filtering.rs`) pins all six values
  from a policy whose two booleans differ, so a source transposition cannot pass;
  `provenance_parameters_spell_an_off_filter_rather_than_omitting_it` pins the `"off"`
  spellings and the floor's conditional absence; the `u64::MAX` length's digits-as-string
  fallback has its own assertion.
- The readable-body test pins `observation-reach-ceiling-bp = 4000`, the `[[read-group]]` row
  keys and first-row values, and all six read-filter keys; the widest-number test now also
  holds the reach ceiling at exactly `i64::MAX`.

## Validation

In the dev container (2026-09-03): `cargo fmt --check` clean;
`cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --all-targets --all-features` — all 16 binaries pass, lib suite **6,055 passed,
0 failed** after the review fixes (`ng::psp` alone: 417; `ng::read::filtering`: 21).

## Assumptions recorded

- The ceiling's typed home is `Bp` (u64) though the generator's cap is `u32` — the header
  records a width in bases and does not import the generator's integer width; the rule set
  bounds it at TOML's ceiling.
- `record_read_filters` overwrites its six keys if called twice; last write wins, matching
  `BTreeMap::insert`. The gatherer calls it once (plan B1).
