# ng step 2 — read preparation: types & interfaces

*Status: architecture draft (2026-07-26), companion to the spec
[`../spec/read_preparation.md`](../spec/read_preparation.md) (the design and its rationale) and to the
shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary + step traits) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). `ng_step_interfaces.md` §Step 2 sketches
the trait only; **this doc is the real interface.** Read preparation is a **generic-path-only** step —
the STR path has none (spec §1). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns for types,
verbs for functions; **STR** in prose ↔ `ssr` in code (this step is generic, so `ssr` does not appear).
Signatures are illustrative; the **contract** is the deliverable. Every "why" lives in the spec.*

## Module home

Two files in the existing `read/` module (`module_layout.md` principle 1b — steps 1 and 2 share it):

- **`src/ng/read/mod.rs`** — the `ReadPreparer` trait and `ReadPrepError`. The trait sits in `mod.rs`
  because it is the module's shared contract; the impls sit beside it.
- **`src/ng/read/left_align.rs`** — `LeftAlignPreparer`, the v1 impl, with its `#[cfg(test)]` block.

A folder is not warranted yet, but this step **does** have a bake-off (unlike step 1): further impls —
a pass-through-only null arm, and the gated re-align mode — land as sibling files in `read/`.

*Correction owed to `module_layout.md`:* its tree still names this file `left_align_baq.rs`. BAQ is
deferred sine die (spec §10), so the file is `left_align.rs`, which is what that doc's own "Naming to
confirm" section already says.

## 1. Types

### 1.1 Nothing new to mint

The input is production's `MappedRead` and the output is production's `PreparedRead`, both reused
as-is (spec §3) — so this step seeds **no scalar newtypes** and defines no read type of its own. One
conversion crosses the seam: `MappedRead.ref_id: usize` → `PreparedRead.chrom_id: u32`. It is the same
narrowing step 1 already does, with the same treatment — a value that did not fit would be a corrupt
record, so it fails loudly rather than truncating silently.

### 1.2 Step-2-local types

```rust
/// A fatal, run-level failure. In v1 there is exactly one source: the reference fetch that
/// left-alignment needs. It is NEVER a per-read verdict — a read that yields nothing is
/// `Ok(None)` (spec §7).
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ReadPrepError {
    /// The reference window for an indel-bearing read could not be fetched.
    #[error("reference access failed during read preparation")]
    Reference(#[source] RefSeqError),
}

/// Reused buffers for `LeftAlignPreparer` — one per worker, never per read.
/// Buffers only: nothing here may change a result.
#[derive(Debug, Default)]
pub struct LeftAlignScratch {
    /// The reference window, refilled per indel-bearing read. Reads with no indel never
    /// touch it (spec §5).
    ref_buf: Vec<u8>,
}
```

**No counts type in v1.** Spec §7 wants a per-sample tally of declines by reason; v1 has no decline
reason at all, so a `ReadPrepCounts` now would be an all-zero struct with no variants to count. The
driver counts what it sees. See §6 — this is the one place the settled signature will have to move.

## 2. The interfaces

The trait as it appears in code. The three signature decisions behind it (by value, `Result` around
`Option`, no reference argument) are argued in spec §6 and recorded in §4 below.

```rust
/// Canonicalise the line-up a mapper gave one read, against the reference around that read's
/// own span. Per read, independent of every other read, and independent of every locus.
pub trait ReadPreparer {
    /// Reused buffers — allocated once per worker, never per read. `()` for an impl that
    /// needs none.
    type Scratch: Default;

    /// Prepare one filtered read.
    ///
    /// `Ok(Some(read))` — the prepared read. `Ok(None)` — no usable observation (never in v1).
    /// `Err` — the run is broken; never a per-read verdict.
    fn prepare_read(
        &self,
        read: MappedRead,
        scratch: &mut Self::Scratch,
    ) -> Result<Option<PreparedRead>, ReadPrepError>;
}

/// v1: pass-through + left-alignment. Generic over both the reference it fetches from and the
/// normalizer it calls — neither is chosen by taste, and both stay visible at the use site.
pub struct LeftAlignPreparer<R: RefSeq, N: AlignmentNormalizer = DefaultAlignmentNormalizer> {
    reference: R,
    normalizer: N,
}

impl<R: RefSeq, N: AlignmentNormalizer> LeftAlignPreparer<R, N> {
    pub fn new(reference: R, normalizer: N) -> Self;
}
impl<R: RefSeq> LeftAlignPreparer<R, DefaultAlignmentNormalizer> {
    /// The recommended construction: algorithm 1a, the GATK/production left-aligner.
    pub fn with_default_normalizer(reference: R) -> Self;
}

impl<R: RefSeq, N: AlignmentNormalizer> ReadPreparer for LeftAlignPreparer<R, N> {
    type Scratch = LeftAlignScratch;
    fn prepare_read(&self, read: MappedRead, scratch: &mut LeftAlignScratch)
        -> Result<Option<PreparedRead>, ReadPrepError>;
}
```

**Contract.** Pure per call: the output is a function of the read, the reference, and the impl's own
constructor state, with no hidden mutation beyond `Scratch` — so the result never depends on call
order or thread count, which the cohort's byte-identity guarantee rests on. Content-preserving: bases
and qualities are copied through untouched, and only the CIGAR is rewritten; `alignment_start` does not
move. A `Scratch` that allocates per read is a defect, not a slow path. **The reference is consulted
only for reads whose CIGAR carries an indel** — a read without one never fetches, never touches the
scratch, and cannot fail. Any `RefSeqError` is fatal and returned as `Err`, never degraded to a drop.

## 3. The per-read shape

Two details that decide whether this compiles cleanly, both invisible from the signatures.

**The CIGAR round-trip.** `AlignmentNormalizer::normalize` works on an ng `Alignment`, not on a
`Vec<CigarOp>`, so the CIGAR is moved out of the read and back — no allocation, and no clone:

```rust
let mut alignment = Alignment { reference_offset: 0, cigar: std::mem::take(&mut read.cigar) };
self.normalizer.normalize(&mut alignment, &read.seq, &scratch.ref_buf);
read.cigar = alignment.cigar;
```

`reference_offset` is `0` because the fetched window *starts* at the read's first aligned base — the
offset exists for callers holding a larger stretch, and a non-zero value here would double-count.

**The window.** `[read.pos, read.pos + cigar_ref_span(&read.cigar))`, fetched **uppercased** via
`RefSeq::fetch_into` — ng's deliberate divergence from production's raw bytes (spec §6). The fetched
length must be checked rather than trusted: production's left-aligner fails *safe* on a short window by
leaving the CIGAR untouched, which would silently skip normalization under a model that says the fetch
cannot fail.

`cigar_ref_span` and the `PreparedRead` build are existing code (§5); the only new helper is a
`cigar_has_indel(&[CigarOp]) -> bool` predicate — production's equivalent is private to `indel_norm`.

## 4. Design decisions — decided

Distilled from the spec; see it for the reasoning. Open items are marked `OPEN:` in §6.

- **The step is generic-path only — decided.** The STR path re-aligns every spanning read, so
  canonicalizing the mapper's CIGAR there would be work nothing reads (spec §1). No `type Locus`, no
  `type Prepared` — the path-owned associated types are gone with the STR arm.
- **Static dispatch, generic type parameter, never `Box<dyn>` — decided.** Billions of calls; a
  virtual call per read is a cost this path cannot carry (spec §6). The per-read *mode* is a matched
  enum, not a second dispatch mechanism.
- **The read arrives by value — decided.** Its buffers move into the `PreparedRead` instead of being
  cloned; step 1 already yields owned `MappedRead`s, so it costs the caller nothing (spec §6).
- **`Result<Option<_>, _>`, not `Option` — decided.** A per-read decline and a broken run are
  different outcomes and must not share a channel; step 1 reached the same shape (spec §6, §7).
- **The reference is a field of the impls that need one, not an `Option<R>` on one type — decided.**
  "This preparer needs no reference" stays a compile-time fact. Rejected: a single configurable type
  with `Option<R>`, which would need a defined answer for "no reference, indel-bearing read" (spec §6).
- **The normalizer is a type parameter with a default, not a hardcoded call — decided.** The alignment
  module ships three normalizers and names `DefaultAlignmentNormalizer` precisely so which one runs
  stays visible at the use site; burying `StructuredLeftAligner` inside the preparer would hide it and
  make the normalizer screen unrepeatable through this step. No cost: it is a unit struct.
- **The uppercased (canonical) reference view — decided (owner, 2026-07-26).** An alignment must not
  treat a base differently because a masking tool lowercased it. ng therefore needs **one** reference
  view where production carries two (spec §6, §11).
- **Fetch only for indel-bearing reads — decided.** Left-alignment is the reference's only consumer
  here, and production's reason for fetching on every read (`F1` needed the same window) does not carry
  over, because ng moved `F1` to step 1 (spec §5).
- **Runs in the parallel per-read stage, one `Scratch` per worker — decided.** Production moved this
  work off the serial path deliberately; the walker only consumes `PreparedRead`s (spec §6).
- **One preparer, and one reference accessor, per worker — decided (owner, 2026-07-26). Never
  shared.** Two reasons, and the first is decisive: **the accessor is stateful.** `WindowedRefSeq`
  holds a `RefCell<Option<(ContigId, RawChromReader)>>` — a resident window and a reader that slides
  forward through the contig — and is documented `Send` but **not `Sync`**, "per-worker ownership, like
  the production fetchers" ([ref_seq.rs:490-518](../../../../src/ng/ref_seq.rs)). Second: even where an
  accessor *could* be shared, sharing the **cache** would be the wrong trade — mutating it from several
  threads costs synchronisation and waiting, to save an amount of memory that does not matter.
  Duplication is cheap under either impl, for different reasons: `WindowedRefSeq` keeps only the
  stretch it is walking, so N windows is still small; `ResidentRefSeq` holds an **`Arc` to immutable
  contig bytes**, so what duplicates is the bookkeeping, not the genome (production's shape — a
  per-worker `RawContigRefCache` over a cloned `Arc`-backed repository, consulted only on a contig
  change, [baq_stream.rs:328](../../../../src/pileup/per_sample/baq_stream.rs)).
  **Consequence:** `R` needs no `Sync` bound, and the driver builds the preparer inside `map_init`.
- **`&self` + `Scratch` rather than step 1's `&mut self` — decided.** With the preparer per-worker,
  `&mut self` would also work and would let the buffer be a field. `&self` is kept for two reasons: it
  is ng's idiom for "stateless algorithm + caller-owned buffers", shared with the alignment traits this
  step composes with ([alignment/mod.rs](../../../../src/ng/alignment/mod.rs)); and it stops an impl
  accumulating cross-read state **in its own fields**, which is what the spec's per-read independence
  property forbids (spec §6). It is not a proof of statelessness — the reference behind it may carry a
  resident window through interior mutability, as `WindowedRefSeq` does — but that state is a *cache*,
  and eviction is explicitly "a hint, never a fact the answer depends on"
  ([ref_seq.rs:225-229](../../../../src/ng/ref_seq.rs)).
- **`WindowedRefSeq` is the accessor to use here — decided.** Preparation needs only the read's own
  footprint, and reads arrive coordinate-sorted: exactly what the windowed reader is built for. It
  keeps a sub-range buffer and extends it on demand
  ([raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs)), where `ResidentRefSeq` would hold a
  whole contig per worker. The preparer stays generic over `R: RefSeq`, so this is the driver's choice
  to make; it is the one to make.
- **The driver MUST evict, and eviction is its call, not `prepare_read`'s — decided.**
  `EvictableRefSeq::evict_before` takes `&mut self`, which `prepare_read(&self, …)` cannot provide, so
  the worker that owns the reference evicts it. **This is not optional housekeeping: without it the
  window grows forward to the whole contig and `WindowedRefSeq` becomes `ResidentRefSeq` with extra
  steps.** Since reads are sorted, the call is trivial — `evict_before(read.pos)` after each read is
  what production does in its own hot loop
  ([baq_stream.rs:335-339](../../../../src/pileup/per_sample/baq_stream.rs)).

  **Two things a reader will get backwards.** *Frequent eviction is the cheap regime, not the
  expensive one:* `evict_before` is a `drain` from the front, so it costs O(bytes **remaining**), not
  O(bytes dropped) ([raw_chrom_reader.rs:317](../../../../src/ng/raw_chrom_reader.rs)) — evicting every
  read keeps the window at about one footprint, and evicting rarely makes each call move more. And
  *eviction never re-reads*; it only drops. But **a fetch behind the evicted point does re-read from
  the file** (`prepend_backward`, [raw_chrom_reader.rs:292-295](../../../../src/ng/raw_chrom_reader.rs)),
  so each worker must be fed **monotonically increasing** packets. A packet handed to a worker behind
  its evicted point makes every read in it re-read from disk.

  `OPEN:` the window bounds its memory only if the caller evicts, and `drain` keeps the allocation, so
  a worker never gives back its high-water mark. Whether the reader should bound itself instead —
  chunked internal storage with whole-chunk eviction, or a dead-prefix offset compacted lazily — is a
  `ref_seq` question, not a step-2 one. Measure on peak RSS, not eviction time.
- **No `bench/` this step yet — decided.** The impls that would populate a frontier (the null arm, the
  re-align mode) do not exist; v1's measurement is the parity fixture (§7).

## 5. Reconciliation with existing code

Every row read before it was written. "reuse" means call it; "port" means ng owns a copy.

| ng name | existing code | action |
|---|---|---|
| `PreparedRead` (output) | [pileup/walker/mod.rs:236](../../../../src/pileup/walker/mod.rs) | **reuse as-is** — production's raw `u32`/`u8` fields, `qname: Arc<str>` included. Misplaced inside `pileup/walker/`; recorded, not acted on (production frozen) |
| `MappedRead` (input) | [bam/alignment_input.rs:78](../../../../src/bam/alignment_input.rs) | reuse as-is — step 1's output |
| `CigarOp` | [pileup/walker/mod.rs:43](../../../../src/pileup/walker/mod.rs) | reuse as-is (already `pub`; the `alignment/` module reuses it too) |
| building the `PreparedRead` | `prepare_passthrough` [baq_engine.rs:405](../../../../src/pileup/per_sample/baq_engine.rs) | **`pub`, and it does the whole build** — `alignment_end`, `mate_role`, `qname`, `mq_log_err`, adaptor boundary, raw-`qual`→`bq_baq`. Call-vs-port is `OPEN:` (§6) |
| ↳ its inner wiring | `mapped_to_prepared` [baq_engine.rs:410](../../../../src/pileup/per_sample/baq_engine.rs) | **private** — unreachable from ng except through `prepare_passthrough` |
| the normalizer trait | `AlignmentNormalizer` [ng/alignment/mod.rs:618](../../../../src/ng/alignment/mod.rs) | reuse — **built**; no `Scratch` on it by design |
| the default normalizer | `DefaultAlignmentNormalizer` = `StructuredLeftAligner` [ng/alignment/mod.rs:655](../../../../src/ng/alignment/mod.rs), [left_align_structured.rs:65](../../../../src/ng/alignment/left_align_structured.rs) | reuse — **built**; a unit struct wrapping production's `left_align_indels` |
| the normalizer's input type | `Alignment` [ng/alignment/mod.rs:77](../../../../src/ng/alignment/mod.rs) | reuse — `{ reference_offset, cigar }`; §3's round-trip |
| the reference | `RefSeq::fetch_into` [ng/ref_seq.rs:146](../../../../src/ng/ref_seq.rs) | reuse — the **canonical** view, not `RawRefSeq` (§4) |
| reference impls | `ResidentRefSeq` [ng/ref_seq.rs:358](../../../../src/ng/ref_seq.rs), `WindowedRefSeq` [ng/ref_seq.rs:500](../../../../src/ng/ref_seq.rs) | reuse — **one instance per worker**, never shared: `WindowedRefSeq` carries a sliding reader in a `RefCell` and is `Send`, not `Sync` (§4) |
| eviction | `EvictableRefSeq::evict_before` [ng/ref_seq.rs:230](../../../../src/ng/ref_seq.rs) | not used in v1; `&mut self`, so driver-level (§4) |
| `ReadPrepError::Reference` | `RefSeqError` [ng/ref_seq.rs:39](../../../../src/ng/ref_seq.rs) | wrap — `#[non_exhaustive] thiserror`, same shape as `ReadFilterError::Reference` [ng/read/filtering.rs:583](../../../../src/ng/read/filtering.rs) |
| reference span of a CIGAR | `cigar_ref_span` [bam/alignment_input.rs:947](../../../../src/bam/alignment_input.rs) | call directly (`pub(crate)`, same crate) |
| `cigar_has_indel` | `is_indel` [indel_norm.rs:79](../../../../src/pileup/walker/indel_norm.rs) | **new** — production's is private; a two-line ng predicate |
| the per-read fold (model) | `process_read` [read_processor.rs:165](../../../../src/pileup/per_sample/read_processor.rs) | model only — ng takes its `F3` stage; `G2`/`F1` are step-1 filters, BAQ is deferred |
| the outcome enum (model) | `ReadOutcome` / `DropReason` [read_processor.rs:129-146](../../../../src/pileup/per_sample/read_processor.rs) | model for the decline-reason shape v1 does not yet need (§6) |
| the parallel deployment (model) | the rayon `map_init` in [baq_stream.rs:331](../../../../src/pileup/per_sample/baq_stream.rs) | model — per-worker reference cache + engine, one packet per contig |

## 6. Open items

- `OPEN:` **`Ok(None)` carries no reason, but spec §7 wants a tally by reason.** Harmless in v1, where
  nothing declines. The first decline reason (BAQ's skip, or a re-align that cannot place a read)
  forces `Option<PreparedRead>` to become a two-variant outcome carrying the reason — production's own
  `ReadOutcome { Prepared, Dropped(DropReason) }` is that shape. Settle it when the first decline
  arrives, not before; it is a signature change with no caller today.
- `OPEN:` **Call `prepare_passthrough` or port it.** Leaning: call it in the first slice, port if its
  avoidable `read.qual.clone()` shows up in a profile (spec §11).
- **Which accessor a worker gets — `WindowedRefSeq` is the fit, but the choice is the driver's.**
  Preparation needs only the read's own footprint and the reads arrive coordinate-sorted, which is
  precisely what the windowed reader is for: it keeps a sub-range buffer, extends it on demand, and
  drops what the walk has passed on `evict_before`
  ([raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs)). `ResidentRefSeq` holds a whole
  contig instead — cheaper per fetch, far more memory. The preparer is generic over `R: RefSeq` and
  takes either, so this is a driver-level memory decision, not a step-2 one (§4).
- **Impl-time confirmation — `Bp` vs raw `u64` at the fetch boundary.** `cigar_ref_span` returns `u32`
  and `RefSeq::fetch_into` takes `u64`; the widening is infallible. Pin the spelling when coding.

## 7. Test & bench shape

- **Unit tests beside the code** (`read/left_align.rs`): a read with no indel round-trips unchanged and
  performs no fetch; a read with an indel in a homopolymer comes back left-aligned; `alignment_start`
  never moves; a fetch failure surfaces as `Err`, not `Ok(None)`.
- **The parity fixture — the regression anchor** (spec §9). Compare against production's
  `process_read(read, None, &mut raw_ref, &cfg)` in-process with `max_read_mismatch_fraction: None`
  (ng moved `F1` to step 1; leaving it on diverges the keep-sets for unrelated reasons). It proves the
  window fetch and the field wiring — **not** left-alignment, which is already byte-parity-checked
  where the normalizer landed ([alignment_normalization.md](../impl_plan/alignment_normalization.md)
  step B1). Follow the `#[cfg(test)]` parity-module pattern of
  [delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs) so shipping ng code keeps no
  test-only dependency on production.
- **Two fixtures, deliberately.** On an **all-uppercase** reference the output must be byte-identical
  to production. On a **soft-masked** reference it must *not* be, and the difference is the size of the
  production defect ng is fixing (spec §6, §9). A single fixture cannot show both.
- **No `bench/`** — see §4.
