# ng generic locus generator — the port (plan 2 of 3)

**Status:** draft, 2026-07-28. Copy production's pileup walker into ng and **prove it computes
exactly what production computes** — before a line of it is changed. Design settled in
[`locus_generation_pileup.md`](../spec/locus_generation_pileup.md) §3 (spec) and
[`../arch/locus_generation_pileup.md`](../arch/locus_generation_pileup.md) *Module home* (types &
interfaces). This turns that design into build order; it is **not** a place for new design.

Follows [`locus_generation_pileup_prerequisites.md`](locus_generation_pileup_prerequisites.md);
precedes [`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md).

---

## Scope

**In:** `src/ng/read/prepared_read.rs` (ng's own read type); the seven copied walker files under
`src/ng/locus_generation/pileup/`; the `RefSeq` → `MultiChromRefFetcher` shim; and the **stage-1
differential** that is this plan's whole point.

**Out (plan 3):** every behaviour change — the no-fill builder, `read_coverage`, REF-only widening,
the read-group split, the per-record counters, the generator wrapper, the region walk, stage-2
parity, the dump tool. **Nothing in this plan changes what the walker computes.**

**Out (never, here):** any edit to `src/pileup/` — production is frozen and this plan copies rather
than reaches in (spec §3).

## Principles (how the order was chosen)

- **Transcribe first, change second.** The rule that paid three times on this branch (`delimit_read`
  over 200,000 randomized cases, `left_align_indels`, the SplitMix64/reservoir port). A copy that is
  *provably* production is the baseline every later change is measured against; without it, plan 3's
  divergences cannot be told from transcription slips.
- **Types first, then implementation** — the read type (A1) before the code that names it (A2).
- **Verify against ground truth.** The oracle is production's own walker over the same input, not
  self-consistency.
- **A differential that passes immediately is suspect.** Milestone B does not end when the
  differential is green; it ends when the fixture has been **shown to discriminate**, by mutating
  each behaviour in turn and watching it fail. This is the project's recurring lesson — a test that
  cannot fail is the failure mode, not the exception.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at
  completion.

## Preconditions (already in place)

- **Plan 1 is complete** — the region stream is owned (its Milestone A) and the shared locus type is
  final (its Milestone B). The type matters here only because plan 3 fills it; the stream matters
  not at all yet.
- **Production's walker is intact and frozen**: `driver.rs`, `open_record.rs`, `cigar_cursor.rs`,
  `decompose.rs`, `active_read_set.rs`, `chain_id_allocator.rs`, `errors.rs` (~5,495 lines).
- **Its test suite exists** — **44** end-to-end tests in [`walker/tests.rs`](../../../../src/pileup/walker/tests.rs),
  plus **69** inline across the seven files (70 `#[test]` markers; `subtract_contribution`.s debug and
  release pair is mutually exclusive by `cfg`), so **113 inherited tests** in any one profile. *(Both
  documents said "46" until 2026-07-29; the number was counted during A4 and corrected here.)*
  `MockFasta`, `snp_read` and `paired_snp_reads` are reusable from ng.s tests (`pub(crate)` under
  `#[cfg(test)]`). Spec §12 classifies them.
- **`CigarOp`, `PileupRecord`, `AlleleObservation`, `AlleleSupportStats` are `pub`** and unchanged by
  ng — reused, not copied.

---

## Milestone A — the copy

- ✅ **A1 — ng's `PreparedRead`.** `src/ng/read/prepared_read.rs`: `PreparedRead`, `MateRole`,
  `ReadLengthError`, copied from `pileup/walker/mod.rs` and extended with
  `read_group: ReadGroupId`. `ReadPreparer` returns this type; `LeftAlignPreparer` threads the group
  through from `AlignedRead`. Keep it **not** `#[non_exhaustive]`, for production's stated reason —
  a new field should break every construction site. *Depends:* —. *Source:* spec §6, arch *Module
  home*.
- ✅ **A2 — the seven files, verbatim.** Copy into `src/ng/locus_generation/pileup/`, renaming
  `driver.rs` → `genome_walk.rs` (the only one of the seven named for a role rather than for what it
  owns). The **only** edits permitted: the module paths, and `PreparedRead` resolving to ng's. Still
  emits `PileupRecord`. Their inline `#[cfg(test)]` modules come along — including `decompose`'s
  oracle, which is what the cursor is parity-tested against. *Depends:* A1. *Source:* spec §3, arch
  *Module home*.
- ✅ **A3 — the reference shim.** `RefSeqFetcher<R: RefSeq>` implementing `MultiChromRefFetcher`.
  Semantically empty: both contracts are canonical uppercase `{A,C,G,T,N}`, verified in the
  implementation and not just the doc. *Depends:* A2. *Source:* arch §1.3.
- ✅ **A4 — the copied suite is green.** All **113** inherited tests (44 end-to-end + 69 inline) pass **unmodified**. Anything needing
  a touch here is a transcription error, not a design change — spec §12 is explicit that this is the
  gate. *Depends:* A3. *Source:* spec §12.

> **Checkpoint A: the copy compiles, and production's own tests pass against it untouched.**
> Pause for review.

---

## Milestone B — the differential, and proving it can fail

- ✅ **B1 — the stage-1 harness.** `pileup/parity.rs`, `#[cfg(test)]`, in the `delimit_parity` /
  `left_align_parity` shape: build one `Vec<PreparedRead>`, drive
  `crate::pileup::walker::run` and ng's `run` with the same `WalkerConfig` and the same fetcher, and
  assert the two `Result<PileupRecord, WalkerError>` streams are equal element for element, plus
  `RunSummary`. Byte-identity is well defined because `PileupRecord`'s `PartialEq` compares the two
  `f32`s **by bits**, so the `NaN` placeholders compare equal. Feed **one** prepared stream to both —
  preparing separately would inject step 2's uppercase divergence. *Depends:* A4. *Source:* spec §3.
- ✅ **B2 — prove the harness discriminates. Own commit, do not bundle.** Mutate each of five
  behaviours in the ng copy in turn and require the differential to **fail**: mate-overlap
  reconciliation, adaptor masking, record widening, the subtract-then-add re-fold, and the column
  depth cap. A fixture that cannot fail is worth nothing, and this is where that is established —
  the fixture is fixed here and reused for the rest of the port. *Depends:* B1. *Source:* spec §12,
  §13.1.
- ✅ **B3 — at scale.** Run the differential on GIAB HG002 and a tomato CRAM under the
  `PVC_PARITY_CASES` convention. Zero divergences, or the port is not done. *Depends:* B2.
  *Source:* spec §13.1.

> **Checkpoint B: ng's walker is production's walker, demonstrably — and the demonstration has been
> shown capable of failing.** Pause for review. **This is the baseline plan 3 measures against; it
> cannot be reconstructed later**, because plan 3's first commit makes the two walkers differ on
> purpose.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | production.s **113** inherited walker tests (44 end-to-end + 69 inline), **unmodified**, green against the copy |
| B1 | element-wise equality of two `PileupRecord` streams + `RunSummary`, one input stream, bitwise `f32` comparison |
| B2 | five mutations, five failures — the fixture is shown to discriminate before it is trusted |
| B3 | zero divergences on GIAB HG002 and a tomato CRAM at `PVC_PARITY_CASES` scale |

## Out of scope (next plan)

All of [`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md): the
no-fabrication rule, `read_coverage`, REF-only widening, the read-group split, the per-record
counters, `PileupGenerator`, the region walk with its halo and clamp, the allocator reset, stage-2
projection parity, the dump tool and its six new fixtures, and the throughput measurement.

**One thing to carry forward rather than rediscover:** the stage-1 differential **dies** when plan 3
lands, because the two walkers then differ by design. What survives as a permanent regression anchor
is narrower — loci where every folded read witnessed the whole footprint must agree with production
forever (spec §3). Plan 3 builds that; this plan builds the thing it is derived from.
