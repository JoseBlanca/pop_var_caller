# ng read preparation (step 2) — implementation plan

**Status:** draft, 2026-07-26. The build order for **step 2, read preparation**: the
`ReadPreparer` trait, its error, and `LeftAlignPreparer` — pass-through + left-alignment, producing
production's `PreparedRead`. Design is settled in [`../spec/read_preparation.md`](../spec/read_preparation.md)
(the what and why) and [`../arch/read_preparation.md`](../arch/read_preparation.md) (types &
interfaces), under the shared arch docs ([step interfaces](../arch/ng_step_interfaces.md),
[module layout](../arch/module_layout.md)). This turns that design into build order; it is **not** a
place for new design — the one item still open (call vs port `prepare_passthrough`) carries a leaning
in spec §11 and is followed here, not re-decided.

This is the plan [`read_filtering.md`](read_filtering.md) deferred ("step 2, read preparation — the
`read/` sibling; its own plan").

---

## Scope

**In:** `src/ng/read/mod.rs` gains the `ReadPreparer` trait and `ReadPrepError`; a new
`src/ng/read/left_align.rs` holds `LeftAlignScratch` and `LeftAlignPreparer<R, N>` — the conditional
reference fetch, the CIGAR round-trip through `AlignmentNormalizer`, and the `PreparedRead` build; a
`cigar_has_indel` predicate; and the parity fixtures, on an uppercase **and** a soft-masked reference.

**Out (later plans, nothing dropped):**

- **The re-align mode (algorithm 2)** — gated on both an unbuilt affine aligner
  ([`alignment_best_path.md`](alignment_best_path.md) Milestone E, itself gated on this step's §4
  trigger question) and an undecided trigger. Spec §4, §10.
- **BAQ** — deferred sine die; when it returns it is a config mode on `LeftAlignPreparer`, not a
  second impl. Spec §10.
- **The pass-through-only null arm** — a five-line impl whose only purpose is a measurement ("does
  left-alignment change calling at all"), which needs a pileup to run it through. It lands with the
  first driver, not here. Arch §2.
- **The driver: who calls `prepare_read`, and the eviction discipline it owes** — the parallel
  per-read stage and its `evict_before` calls belong to the ng pileup, which does not exist. Arch §4.
- **A `ReadPrepCounts` tally** — v1 has no decline reason to count. It arrives with the first one
  (BAQ's skip, or a re-align that cannot place a read), together with the outcome-enum change arch §6
  records. Spec §7.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The cheap path before the expensive one.** The no-indel path (Milestone B1) touches no reference
  and cannot fail; it is built and tested first, so the fetch/normalize path lands against a working
  build rather than alongside one.
- **Isolate the step whose failure is silent.** B2 is the only step here that can produce a *quietly
  wrong* answer — a mis-sized window or a non-zero `reference_offset` yields a mis-placed indel, which
  is a wrong variant, not a panic. It lands as its own commit with its oracle green before and after.
- **Reuse over rewrite.** The shifting is production's `left_align_indels`, already wrapped and
  byte-parity-checked as `StructuredLeftAligner`; the `PreparedRead` field wiring is production's
  `prepare_passthrough`. This step supplies the window, the round-trip, and the fetch policy — nothing
  else, and no alignment logic is re-derived.
- **Verify against ground truth.** The north-star test is byte-parity with production's `process_read`
  under `--no-baq` on an **uppercase** reference — and a *deliberate* divergence on a soft-masked one,
  which is the measurement, not a failure (spec §9).
- **Incremental, with pauses.** One milestone, then stop for review.
- **Ungated / container builds.** `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place — verify before A1)

- **The normalizers are built.** [`alignment_normalization.md`](alignment_normalization.md) A1–D1 are
  ✅: `AlignmentNormalizer`, `StructuredLeftAligner` (1a, byte-parity vs production's
  `left_align_indels`), and `DefaultAlignmentNormalizer = StructuredLeftAligner`
  ([src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs)).
- **The reference access is built** — `RefSeq::fetch_into` (the canonical view this step uses),
  `InMemoryRefSeq` (the test oracle), `ResidentRefSeq`, and `WindowedRefSeq` with `evict_before`, all
  in [src/ng/ref_seq.rs](../../../../src/ng/ref_seq.rs). **Note:** [`foundations.md`](foundations.md)
  Milestone C is still ☐ in that plan, but the code exists and the typed-region walk already calls
  `evict_before` on it — the checkboxes are stale, not the dependency. Confirm by reading, then fix
  those boxes.
- **Step 1 is done** and yields **owned** `MappedRead`s (`Item = Result<MappedRead, ReadFilterError>`),
  which is what lets `prepare_read` take the read by value at no cost.
- **The production reuse targets and the oracle** — `PreparedRead`
  ([pileup/walker/mod.rs:236](../../../../src/pileup/walker/mod.rs)), `prepare_passthrough` (`pub`,
  [baq_engine.rs:405](../../../../src/pileup/per_sample/baq_engine.rs)), `cigar_ref_span`
  (`pub(crate)`, [bam/alignment_input.rs:947](../../../../src/bam/alignment_input.rs)), and
  `process_read` ([read_processor.rs:165](../../../../src/pileup/per_sample/read_processor.rs)) — the
  parity oracle.
- **Fixture builders** in [read/input/test_fixtures.rs](../../../../src/ng/read/input/test_fixtures.rs)
  (BAM writer + fixture FASTA) — Milestone C adds a lowercased copy of that FASTA, nothing more.

---

## The steps

### Milestone A — the trait, the error, the impl types (types, no logic)

**✅ A1. `ReadPreparer` + `ReadPrepError` in `read/mod.rs`.**
The trait (`type Scratch: Default`; `prepare_read(&self, read: MappedRead, scratch) -> Result<Option<PreparedRead>, ReadPrepError>`)
and the `#[non_exhaustive] thiserror` error with its single `Reference(RefSeqError)` variant; declare
`pub mod left_align;`. Test: a `#[cfg(test)]` stand-in impl driven **through a generic bound**, which
is the only place `Scratch: Default` and the static-dispatch rule actually bite — a concrete call
compiles either way. *Source:* arch §1.2, §2; spec §6, §7.

**☐ A2. `LeftAlignPreparer` + `LeftAlignScratch` types.**
The struct (`reference: R`, `normalizer: N`), `new`, `with_default_normalizer`, and the scratch
holding `ref_buf: Vec<u8>`. No trait impl yet. *Depends:* A1. *Source:* arch §2.

> **Checkpoint A:** the trait compiles, is implementable through a generic bound, and the error
> surface is fixed. Pause for review.

### Milestone B — the transform (the heart)

**☐ B1. `cigar_has_indel` + the no-indel path.**
The predicate (production's `is_indel` is private to `indel_norm`, so this is ng's two-liner), then
`impl ReadPreparer for LeftAlignPreparer`: a read whose CIGAR carries no indel is built straight from
the record — **no fetch, no scratch touched, cannot fail**. Build via `prepare_passthrough` (spec
§11's leaning; a port is the fallback if its `read.qual.clone()` ever profiles). Tests: the read
round-trips with its CIGAR and qualities intact and every `PreparedRead` field wired; and the
reference is *never consulted* — assert it with a reference impl that fails the test if fetched.
*Depends:* A2. *Source:* spec §5, arch §2, §3.

**☐ B2. The indel path. Own commit; do not bundle.**
Fetch `[read.pos, read.pos + cigar_ref_span(cigar))` **uppercased** via `RefSeq::fetch_into` into
`scratch.ref_buf`; **check the returned length rather than trusting it**; wrap the CIGAR in
`Alignment { reference_offset: 0, cigar: mem::take(&mut read.cigar) }`, normalize, move it back, build.
Tests: an indel in a homopolymer comes back left-aligned; `alignment_start` never moves; a fetch
failure surfaces as `Err`, never `Ok(None)`. **Why isolated:** every failure mode here is silent — a
mis-sized window or a non-zero `reference_offset` mis-places an indel, which is a wrong variant with no
crash and nothing in the output to show for it. Its oracle is production's `left_align_indels` on the
same read and window (the *shifting* is already byte-parity-checked, so what this step risks is the
**window**). *Depends:* B1. *Source:* spec §5, §6; arch §3.

> **Checkpoint B:** both paths work against `InMemoryRefSeq`; the no-indel path provably never
> fetches; an indel shifts and the placement start holds. Pause for review.

### Milestone C — parity, and the measurement

**☐ C1. The uppercase parity fixture — the port anchor.**
In-process against `process_read(read, None /* no BAQ */, &mut raw_ref, &cfg)` with
**`max_read_mismatch_fraction: None`** (ng moved `F1` to step 1; leaving it on diverges the keep-sets
for reasons unrelated to preparation). On an **all-uppercase** reference every field of every
`PreparedRead` must be byte-identical. Follow the `#[cfg(test)]` parity-module pattern of
[delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs), so shipping ng code keeps no
test-only dependency on production. *Depends:* B2. *Source:* spec §9, arch §7.

**☐ C2. The soft-masked fixture — the divergence, measured.**
The same BAM against a **lowercased** copy of the fixture FASTA. Production's left-alignment stalls
there (it compares raw bytes against uppercased read bases); ng's does not. The test asserts the
divergence **exists** and reports how many reads' CIGARs differ — that count is the size of the
production defect ng is fixing, and it is the first number attached to a claim the spec marks "not yet
measured". *Depends:* C1. *Source:* spec §6, §9, §11.

**☐ C3. The `F1`-invariance assertion.**
Across the C1 fixture, assert `read_exceeds_mismatch_fraction` returns the same verdict before and
after left-alignment. This is the check spec §5 names as what would settle its step-1/step-2 ordering
argument, which currently rests on reasoning rather than a debug assertion that compiles out.
*Depends:* C1. *Source:* spec §5.

> **Checkpoint C:** byte-parity on uppercase, a measured divergence on soft-masked, and the `F1`
> ordering argument discharged. **Step 2 is complete.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the trait is implementable **through a generic bound** (where `Scratch: Default` and static dispatch bite); the error surface compiles |
| B | unit tests against `InMemoryRefSeq`: no-indel round-trip with the reference provably untouched; an indel left-aligned with `alignment_start` fixed; a fetch failure as `Err`, not `Ok(None)` |
| C | **byte-parity vs production's `process_read` (`--no-baq`, `F1` off) on an uppercase reference**; a *required* divergence on a soft-masked one, with the differing-read count reported; `F1` verdict invariant across left-alignment |

## Out of scope (next plans)

- **The re-align mode** — [`alignment_best_path.md`](alignment_best_path.md) Milestone E (the affine
  aligner) plus the trigger question spec §4 leaves open. That milestone's gate points back at this
  step, so if it is still open when the aligner is wanted, the milestone moves into a follow-on plan
  here (spec §11).
- **BAQ** — a later config mode on `LeftAlignPreparer`, reusing `BaqEngine::process` (spec §10).
- **The pass-through-only null arm and its measurement** — with the first driver.
- **The ng pileup / driver** — the parallel per-read stage, per-worker preparer construction, and the
  `evict_before` discipline a windowed reference needs (arch §4).
- **`WindowedRefSeq`'s self-bounded memory** — whether the window should bound itself (chunked storage
  with whole-chunk eviction, or lazy compaction) instead of trusting the caller to evict. A `ref_seq`
  question, measured on peak RSS (arch §6).
