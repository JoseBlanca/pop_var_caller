# ng generic locus generator — prerequisites (plan 1 of 3)

**Status:** draft, 2026-07-28. Two enabling changes to **already-built, already-reviewed** modules
that the generic locus generator needs before its first line is written. Design settled in
[`locus_generation_pileup.md`](../spec/locus_generation_pileup.md) (spec) and
[`../arch/locus_generation_pileup.md`](../arch/locus_generation_pileup.md) (types & interfaces).
This turns that design into build order; it is **not** a place for new design.

Plans 2 and 3: [`locus_generation_pileup_port.md`](locus_generation_pileup_port.md) (copy the
walker, prove it identical) and
[`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md) (change it, wrap it,
measure it).

---

## Scope

**In:** the two changes to built modules, each on the `bundle_threshold`-rename model — its own
cargo-verified commit, its own review, never a drive-by inside the generator's work.

1. **`src/ng/read/input/` returns an owned region stream** — three borrows become `Arc`s so a
   generator can hold the stream across `next_locus` calls (spec §2, arch §2.2).
2. **The four shared locus-type changes** — `ReadCoverage`, and three fields on `ObservedSequence`
   (spec §10). They land together because they share **one** STR fixture rebaseline.

**Out (later plans):**

- Everything in `src/ng/locus_generation/pileup/` — plans 2 and 3.
- Moving `PreparedRead`/`CigarOp` out of `pileup/walker/` — the recorded misplacement, deferred to
  the port-back (spec §10).
- Any change to `src/pileup/`, `src/psp/`, `src/var_calling/`, `src/vcf/` — production is frozen and
  this work touches none of it (spec §3).

## Principles (how the order was chosen)

- **Representational first, semantic second.** A1–A2 and B1 change *how* things are spelled and must
  move no data; B2 changes *what is emitted* and rebaselines fixtures. Landing them in that order
  means a fixture diff at B2 is unambiguously B2's.
- **One rebaseline, not four.** All four shared-type changes touch types both generators fill, so
  each done separately rebaselines the STR fixtures again (spec §10).
- **Verify against ground truth.** A's oracle is the existing BAM/CRAM parity test — the change is
  representational, so **no read may move**. B1's is the STR dump before/after, equivalent modulo
  the coverage encoding.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at
  completion.

## Preconditions (already in place)

- **`main` is merged** — this branch is at `eb2857c`, carrying the read-fetch perf work and the
  read-group work (`AlignedRead`, `read_groups.rs`). Both are load-bearing: A builds on the borrow
  shape that merge left, B2 on `ReadGroupId`.
- **The read-input suite and the BAM/CRAM parity oracle exist** and are green — they are A's oracle.
- **The STR generator and its fixtures exist** ([`ssr.rs`](../../../../src/ng/locus_generation/ssr.rs))
  — they are B's oracle and B2's rebaseline target.
- Two standard validation commands are red for unrelated reasons and are excepted by hand — see
  PROJECT_STATUS *Standing project-wide items*.

---

## Milestone A — an owned region stream

The obstacle: `reads_in_region` returns `SampleRegionReads<'_, R>`, and `LocusGenerator` lends
`&SampleReads` per call with no lifetime parameter, so a generator cannot hold the stream between
`next_locus` calls (arch §2.2). Three borrows have to go.

- ✅ **A1 — `Arc<sam::Header>`.** `AlignmentFile` holds its header behind an `Arc` and hands clones
  out; `BamRegionSource` / `CramRegionSource` hold `Arc<sam::Header>` instead of `&'a sam::Header`.
  An independent `Arc`, **not** a reference into an `Arc`'d file — that is what keeps nothing
  self-referential. *Depends:* —. *Source:* arch §2.2.
- ☐ **A2 — `Arc<AlignmentFile>`, and the third borrow.** `SampleReads.files: Vec<Arc<AlignmentFile>>`;
  `BorrowedReader` / `RegionReads` / `RegionSource` hold `Arc<AlignmentFile>`; and
  **`resolution: &'a ReadGroupResolution` becomes owned or `Arc`'d** — it arrived with the
  read-group merge and is the borrow the arch's first draft missed. `reads_in_region` then returns a
  stream with no lifetime. *Depends:* A1. *Source:* arch §2.2, §4 P1.
- ☐ **A3 — prove nothing moved.** The BAM/CRAM parity oracle and the whole read-input suite pass
  unchanged; add one test that a returned stream **outlives** the `&SampleReads` borrow that made
  it, which is the property the generator needs and the only new one. *Depends:* A2. *Source:*
  arch §4 P1.

> **Checkpoint A: the stream is owned and no read moved.** Pause for review.

---

## Milestone B — the shared locus type

Four changes to types both generators fill. B1 is representational; B2 changes STR output.

- ☐ **B1 — `ReadCoverage` → `Complete` + `Observed { offset_in_locus, positions_covered }`.**
  **Own commit, do not bundle** — a coverage encoding that is subtly wrong is a wrong depth, not a
  panic. Covers: the enum; `num_obs_along_locus`'s arm; the STR generator's four minting sites (two
  of which pass the variant **as a function value**, so that helper is restructured, not retyped),
  its complete/partial tally and its sort key; and the six dump tools `--all-targets` builds.
  *Oracle:* the STR dump before and after must be **equivalent modulo the encoding** —
  `PartialLeft(n)` ⇔ `Observed { 0, n }`, `PartialRight(n)` ⇔ `Observed { len - n, n }` — with row
  counts and every support field byte-identical. *Depends:* —. *Source:* spec §6, §10.1.
- ☐ **B2 — `ObservedSequence` gains `read_group: ReadGroupId` and `placed_left: u32`.** The STR
  generator fills both; its rows now split by group, which is the rebaseline. `placed_start` is
  **not** added. *Oracle:* on a single-read-group fixture the STR dump is byte-identical to B1's; on
  a multi-group fixture rows split and their `num_obs` sum to the single-group totals. *Depends:*
  B1. *Source:* spec §6, §10.2, §10.3.
- ☐ **B3 — the doc fold-ins.** [`locus_generation.md`](../spec/locus_generation.md) §3 and
  [`../arch/locus_generation.md`](../arch/locus_generation.md) §1 carry the type: record the new
  `ReadCoverage`, the two new fields, and the `reads_without_observation` caveat (its wording is
  broader than the generic path fills). Also [`read_preparation.md`](../spec/read_preparation.md) §3
  — ng now owns its own `PreparedRead`, reversing its reuse-as-is decision. **Docs only, no code.**
  *Depends:* B2. *Source:* spec §10.

> **Checkpoint B: the shared type is final, and the STR fixtures were rebaselined exactly once.**
> Pause for review — plan 2 cannot start until the type it fills is settled.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the existing BAM/CRAM parity oracle + the read-input suite, **unchanged** — a representational change moves no read; plus one new test that the stream outlives the borrow |
| B1 | the STR dump equivalent modulo the coverage encoding, row counts and support fields byte-identical |
| B2 | single-group STR dump byte-identical to B1; multi-group rows split and sum back to the single-group totals |
| B3 | docs only — reviewed, not tested |

## Out of scope (next plans)

- **Copying the walker and proving it identical** — [`locus_generation_pileup_port.md`](locus_generation_pileup_port.md).
- **The behaviour changes, the generator, the measurements** — [`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md).
- **Whether the STR path should split its rows by group at all** — the STR cohort work's call, handed
  over in spec §10.2; this plan only makes it expressible.
