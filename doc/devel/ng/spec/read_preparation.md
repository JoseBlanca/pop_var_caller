# ng — read preparation (the generic SNP/indel path)

*Status: design spec, rewritten 2026-07-25. **This is the single read-preparation spec** — read
preparation is a **generic-path-only** step (§1 says why), so there is no shared preamble and no STR
sibling; the former `read_preparation_generic.md` and `read_preparation_ssr.md` are short redirects
here. **What is built, and what v1 is:** the left-alignment transform ships (`AlignmentNormalizer` +
three impls in `src/ng/alignment/`); **BAQ is deferred sine die** (§10 — not of interest now); the
**re-align mode's aligner (algorithm 2, affine) and its trigger are both unbuilt/gated** (§4, §10); and
the **`PreparedRead`-producing preparer is not built**. So **v1 read preparation is pass-through +
left-alignment** — nothing else. Grounded in the production `process_read` fold
([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)). The code-facing companion
is [`../arch/read_preparation.md`](../arch/read_preparation.md) (types, trait signatures, the
reconciliation table). Naming: **STR** in prose, `ssr` in code.*

---

## 1. What read preparation is — and why only the generic path has it

Read filtering decided *which* reads survive ([`read_filtering.md`](read_filtering.md)). Read
preparation **canonicalises the line-up the mapper gave a read**: it rewrites a filtered `MappedRead`
against the reference around the read's own span, producing a `PreparedRead` that is still a read.

**Why only the generic path has it: the STR path throws the mapper's line-up away.** It re-aligns
every spanning read against `flank + tract + flank` with the repeat-aware delimiter
([`alignment.md`](alignment.md) §4.2), so canonicalizing the CIGAR the mapper produced would be work
whose result nothing reads. The mapper's CIGAR survives there only as **coordinate arithmetic** — the
read's footprint and the reference-window-to-read-offset mapping used to slice out the bases over the
locus window ([ssr.rs](../../../../src/ng/locus_generation/ssr.rs) `read_footprint` / `ref_to_read`);
the measurement itself comes from the re-alignment, which never consults it. Left-aligning first would
therefore *shift which bases get sliced out* and gain nothing. So **the STR path has no read
preparation** — it goes filtering → observation generation
([`locus_generation_ssr.md`](locus_generation_ssr.md)), with the alignment as a component of the latter.

That is also why the STR per-read operation is not "preparation" even by name: its output is an
**observation about one locus**, not a canonicalised read — the same read at another tract comes out
differently. An earlier design forced both paths under one `ReadPreparer` trait with a path-owned
output and a "the STR path always aligns" mode; it is retired.

**Preparation is locus-independent, and the interface leans on it** — the transform needs no locus, so
one `PreparedRead` serves *every* locus the read overlaps, and `prepare_read` takes no locus argument
(§6). That is a property this step *has*, not the reason the STR path lacks it.

**Non-goals.** Preparation never *drops* a read for a whole-read property — that was filtering. It does
not decompose a read into per-position events, apply the adaptor mask, or reconcile overlapping-mate
qualities — those need the locus-column context and are the **pileup walker's** job (§8, §10). It does
not generate candidate alleles or compute read likelihoods (step 7). It does not recalibrate base
qualities (BQSR). And **it does not reassemble.

In the future we might try re-aligning some reads, but that's not a goal at this moment.

---

## 2. How much work a read needs — the three modes

Preparation picks one of three modes per read. All three produce the **same output type**
(`PreparedRead`, §3), which is what makes them interchangeable and comparable.

| mode | what it does | when | v1? |
|---|---|---|---|
| **pass through** | nothing to the placement | the read shows no insertions or deletions | ✅ |
| **canonicalize** | rewrite the read's indels into their leftmost equivalent spelling | the read has indels and its placement is trusted | ✅ |
| **re-align** | discard the mapper's line-up and compute a fresh one from the read's bases | the read's placement is not trusted (§4) | ❌ gated (§4, §10) |

**Pass-through is a fast path, not a different answer.** Left-alignment shifts indels; a read with no
indels has nothing to shift, so canonicalizing it is provably a no-op. (With BAQ deferred, §10,
pass-through is the whole of "do nothing" — there is no quality-capping step left to decide about.)
**v1 does not implement it as a branch:** the ported normalizer already early-returns on a CIGAR with
no indel, before it touches the reference and without allocating
([indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs) `left_align_indels`;
[left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs) — "the caller does
**not** need to pre-filter no-indel reads"). A mode enum here would only duplicate that scan.

**Canonicalize is about spelling.** The same indel can be written at several equivalent reference
positions when it sits in or near a repeat — the gap slides without changing a base of the result.
Left-alignment picks the leftmost spelling so equivalent variants get an identical one; otherwise the
reads supporting them scatter across several weak candidates instead of pooling into one strong one. The
operation lives in the alignment module ([`alignment.md`](alignment.md) §6) as the `AlignmentNormalizer`
trait (**built** — `StructuredLeftAligner`, the GATK/production port and default; `RepeatedLeftAligner`,
the freebayes-shaped repeated-pass form; `FixpointLeftAligner`, the fail-loud fixpoint wrapper).
Preparation *calls* a normalizer; it does not re-implement left-alignment. **Production also caps base
qualities here (BAQ); ng defers that sine die (§10), so v1 canonicalize is left-alignment only.**

**Re-align is the only mode that questions the mapper — and it is not built.** The other two accept the
read's placement; this one throws it away and computes a new line-up with a general (affine) best-path
alignment algorithm. **Both halves are missing:** the affine aligner is algorithm 2, *gated* and unbuilt
([`alignment.md`](alignment.md) §4.1; `alignment_best_path.md` Milestone E), and its trigger is
undecided (§4). It is the expensive, rare mode and the only route by which a mis-placed read is rescued.

---

## 3. The output — `PreparedRead` (reused from production)

Every mode yields the **existing production `PreparedRead`**
([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)), **reused as-is** — it already *is* the
transform's output, field for field, exactly as step 1 reuses `MappedRead`. Copied verbatim from
production, **not redefined** — the types are production's raw integers, not ng's newtypes, which is
what "reuse" means here:

> **⚠ Fold-in, 2026-07-28 — this decision is REVERSED: ng owns its own `PreparedRead`.** The section
> below still describes production's type accurately and is worth reading as the field inventory,
> but "reuse as-is" no longer holds.
>
> **Why.** A read's library membership is a property of the read, like its MAPQ or its strand, and a
> per-library error model cannot be added later from a tally that merged the libraries together — so
> the prepared read has to carry `read_group: ReadGroupId`. Production's `PreparedRead` has no such
> field, and **production is frozen**, so the group would die at the preparation boundary. An earlier
> draft proposed adding the field to production's type (preceded by relocating `ReadGroupId` to a
> crate-visible home); that is **withdrawn**. ng copies the type into `src/ng/read/` and extends it
> there, which edits production zero times and makes the *next* field free — ng's read type will
> change again (BAQ if it returns, the re-align mode), and each of those would otherwise have been a
> fresh edit to frozen code.
>
> **What it costs and does not cost.** The copy is not an extra cost on top of copying the walker, it
> is what makes copying the walker the right call: every module that looked reusable
> (`cigar_cursor`, `decompose`, `active_read_set`, `chain_id_allocator`) names `PreparedRead` in its
> signatures, so an ng-owned read type reaches all of them whatever else is decided. The preparer
> threads the group through from `AlignedRead`, which already carries it.
>
> *Unchanged by this:* `PreparedRead`'s home inside `pileup/walker/` is a recorded misplacement —
> preparation produces it, the pileup only consumes it — deferred to the port-back. ng copying the
> type neither pays that debt nor worsens it. Source:
> [locus_generation_pileup.md](locus_generation_pileup.md) §3, §6, §10.

```rust
pub struct PreparedRead {
    pub chrom_id: u32,                    // index into the merged ContigList
    pub alignment_start: u32,             // 1-based
    pub alignment_end: u32,               // cached, so the walker never re-walks the CIGAR for span
    pub cigar: Vec<CigarOp>,              // canonicalized (or, later, re-aligned) — unlike MappedRead.cigar
    pub seq: Vec<u8>,                     // uppercase ACGTN
    pub bq_baq: Vec<u8>,                  // in v1: the read's base qualities, passed through UNCAPPED
                                          //   (BAQ deferred, §10). The field is production's; a future
                                          //   BAQ would cap it to min(BQ, BAQ).
    pub mq_log_err: f64,                  // derived ln(P_err) from MAPQ (not on MappedRead)
    pub mapq: u8,                         // raw, preserved
    pub is_reverse_strand: bool,          // decoded from flag
    pub qname: Arc<str>,                  // MappedRead carries Vec<u8>: one Arc<str> allocation per
                                          //   read, with a UTF-8-lossy fallback (qname_to_arc). The
                                          //   walker keys mate pairing on it, so it is not optional.
    pub mate_role: MateRole,              // Solo | FirstOfPair | SecondOfPair (from flag bits)
    pub adaptor_boundary: Option<u32>,    // carried through, applied later by the walker (§10)
}
```

A `PreparedRead` carries **no** per-base overlap adjustment — that is pairwise and happens in the walker
(§10), which is what keeps it self-contained and pairwise-independent. It is deliberately **not**
`#[non_exhaustive]` and has no `Default`, so a caller that misses a field fails to compile rather than
absorbing a wrong value.

Living inside `pileup/walker/` is a misplacement — preparation produces it, the pileup only consumes
it — but production is frozen, so this is **recorded, not acted on**. The same debt is already on
record for `CigarOp`, with the same resolution: porting ng back is the moment to move both
([alignment/mod.rs](../../../../src/ng/alignment/mod.rs) — `Alignment`, "Reusing production's
`CigarOp`").

---

## 4. Choosing the mode — the part that is not settled

Pass-through and canonicalize are chosen from the read's own alignment record: does it carry indels or
not. That needs nothing but the read, and it is all of v1.

**Re-align is different: it needs a judgement the read cannot make about itself.** "The mapper's answer
here is not trustworthy" is a property of the *place*, not one read — a region where reads disagree,
pile up mismatches in one column, or clip at the same offset. Nothing in the current ng step map
produces that judgement. This is one of the two reasons the re-align mode is gated; the
other is that its aligner (algorithm 2) is unbuilt.

---

## 5. The transform detail

Production runs one per-read fold, `process_read`, whose stages are
`G2 bad-CIGAR → F3 left-align → F1 mismatch-fraction → BAQ`. ng assigns the two *rejects* — `G2` and
`F1` — to **step 1** (filters #9 and #8 in `read_filtering.md` §3), and **defers `BAQ` sine die**
(§10). What remains for v1 preparation is `F3`:

- **Indel left-alignment** (the canonicalize mode) — call an `AlignmentNormalizer` (§2). Rewrites
  **only the CIGAR**; bases, qualities and the placement start are untouched. (The default normalizer is
  `StructuredLeftAligner` — the GATK line-up, at one remove: it wraps production's `left_align_indels`,
  which is itself a port of GATK's `AlignmentUtils.leftAlignIndels`
  ([indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs), module note). Production calls it with
  end-deletion stripping *off*, so the alignment's start never moves —
  [left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs), "It does not move
  `reference_offset`".) **Built.** `bq_baq` is filled by copying the read's raw
  qualities through uncapped (§3), production's `--no-baq` behaviour.

**Most reads never touch the reference, and that is a design point, not an optimisation.**
Left-alignment is the *only* consumer of the reference in this step — it needs it intrinsically, since
a gap may slide left only across bases that match **both the reference and the read**
([indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs), module note; the shift core is handed
`[ref_bases, read_seq]`). Everything else in the `PreparedRead` is read off the record: `alignment_end`
from the CIGAR, `mate_role` and strand from the flag, `mq_log_err` from MAPQ, `qname`, the adaptor
boundary. So **fetch only when the read's CIGAR carries an indel**; a read without one is built
straight from the record, with no fetch and no scratch touched.

Production *does* fetch for every read, and the reason does not carry over: `F1` needed the same window
on every read, so it fetched once and shared
([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs) — "F3 + F1 both need the
read's raw reference slice"). **ng moved `F1` to step 1**, which leaves left-alignment alone here.

**Three mechanical details the trait does not show, each a wrong answer if missed:**

- **The reference window is production's exact one** — `[pos, pos + cigar_ref_span(cigar))`, so the
  first reference byte is the read's first aligned base
  ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)). The normalizer takes the
  whole stretch *plus* an offset into it, so with that window the offset is `0`. **ng fetches it
  uppercased** (`RefSeq::fetch_into`), unlike production — see below.
- **The normalizer works on an `Alignment`, not on a CIGAR** — `normalize(&mut Alignment, read,
  reference)`. v1 wraps the read's `cigar` in an `Alignment { reference_offset: 0, cigar }` and takes
  it back out.
- **Check the fetched length; do not trust it.** `left_align_cigar` fails *safe* on a window that does
  not cover the read's footprint: it returns the CIGAR untouched rather than corrupting it
  ([indel_norm.rs](../../../../src/pileup/walker/indel_norm.rs), the under-provisioned-input guard). So
  a short buffer silently skips normalization instead of failing — which under ng's fatal fetch model
  (§7) should never arise, and therefore must be enforced by the fetch rather than left to that guard.

**Why the step-1/step-2 split is safe.** Production runs `F1` *after* left-alignment; ng runs mismatch
filtering in step 1, *before* preparation's left-alignment. The verdict does not move because a gap
only shifts across bases the read and the reference share, so the match/mismatch tally is unchanged —
ng's step-1 code makes the same argument
([filtering.rs](../../../../src/ng/read/filtering.rs) `verdict_post_decode`). **That is an argument,
not a proof, and it is worth marking as one:** production's `left_align_indels` does carry a
mismatch-count invariant, but it is `#[cfg(debug_assertions)]` — out of the build we run — and it
counts raw, case-sensitive mismatches, whereas `F1` counts only mismatches above a BQ floor over
columns where both bases are ATGC
([alignment_input.rs](../../../../src/bam/alignment_input.rs) `read_exceeds_mismatch_fraction`). What
would settle it is one assertion on the parity fixture: the `F1` verdict identical before and after
left-alignment. The bad-CIGAR check (`G2`/#9) sees the raw decoded CIGAR in both.

The re-align transform (algorithm 2) would slot in here as a third per-read outcome when it is built and
triggered (§4); it is out of v1.

---

## 6. The interface — per read, statically dispatched, reused buffers

Three properties every implementation upholds:

- **Per read, and independent of every other read.** No mate, no neighbouring read, no locus-column
  context. This is what makes preparation parallel with deterministic output, and why the genuinely
  pairwise work (reconciling overlapping mates' qualities) sits downstream in the walker. It is also why
  reassembly cannot be a mode here: assembling haplotypes needs every read in a region at once.
- **It re-places what the read already says; it does not invent sequence.**
- **No usable observation is a result, not an error** (§7).

**Dispatch is resolved at compile time.** Preparation runs on every read of every sample — billions of
calls — so which implementation runs is fixed by a generic type parameter, never a trait object
(`Box<dyn …>`): a virtual call per read is a cost this path cannot carry. The per-read *mode* (§2),
which genuinely varies read to read, is a matched enum, not a second dispatch mechanism.

**Buffers are caller-owned and reused.** So preparation threads a reusable scratch value, as the
read-likelihood models do. In `LeftAlignPreparer` it holds **the reference-window buffer**, and that is
the whole reason it exists: `prepare_read` takes `&self`, and `RefSeq::fetch_into` writes into a
caller-owned `Vec<u8>` ([ref_seq.md](ref_seq.md)), so without a scratch every indel-bearing read
allocates one. It is **not** the normalizer's buffers: `AlignmentNormalizer` has no `Scratch` — those
algorithms fill no matrix, and one that wants a buffer owns it as a field
([alignment/mod.rs](../../../../src/ng/alignment/mod.rs)) — and the default normalizer is a unit
struct. The re-align mode's alignment matrices join the scratch when that mode lands.

```rust
pub trait ReadPreparer {
    /// Reused buffers — allocated once per worker, never per read. `()` for an
    /// implementation that needs none.
    type Scratch: Default;

    /// Prepare one filtered read against the reference around its OWN span.
    ///
    /// No locus argument (the generic transform is locus-independent) and **no reference
    /// argument** — an implementation that needs a reference holds its own accessor; one
    /// that does not, does not.
    ///
    /// `Ok(None)` — no usable observation, counted by reason (§7). In v1 nothing returns it.
    /// `Err` — the run is broken (a failed reference fetch). Never a per-read verdict (§7).
    fn prepare_read(
        &self,
        read: MappedRead,
        scratch: &mut Self::Scratch,
    ) -> Result<Option<PreparedRead>, ReadPrepError>;
}
```

**Three decisions in that signature.**

- **The read arrives by value.** Its `seq`, `qual` and `cigar` buffers *move* into the `PreparedRead`
  rather than being cloned per read — production takes it by value for exactly this reason
  ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs), the `--no-baq` arm), and
  step 1 already yields owned `MappedRead`s, so by-value costs the caller nothing.
- **The result is a `Result` around an `Option`,** because the step has two kinds of bad news and §7
  forbids confusing them. A per-read decline is `Ok(None)` and the run continues; a failed reference
  fetch is `Err` and the run is over. An `Option` alone has one channel and could express only one of
  them. Step 1 reached the same shape, for the same reason
  ([read_filtering.md](read_filtering.md) §7).
- **The reference is not in the interface at all** — it is a field of the implementations that need
  one, not an `Option<R>` on a single type. `LeftAlignPreparer<R: RefSeq>` holds `R`; a
  pass-through-only preparer holds nothing and can never fail. This keeps "this preparer needs no
  reference" a compile-time fact rather than a runtime `None` with an undefined case (no reference and
  an indel-bearing read), and it matches both the `Scratch` associated type above and ng's own reason
  for splitting the reference capabilities apart ([ref_seq.md](ref_seq.md)).

The v1 implementation is **`LeftAlignPreparer`** (pass-through + canonicalize). The trait carries **no
`type Locus`** (the generic path needs none) and **no `type Prepared`** (always `PreparedRead`); the
path-owned associated types an earlier sketch had — which could not compile as a `Box<dyn ReadPreparer>`
— are gone with the STR arm. The trait still exists because further impls bench behind it: the re-align
mode, and a pass-through-only preparer as the null arm (does left-alignment change calling at all — a
different question from the normalizer screen, which compared normalizers to each other).

**`LeftAlignPreparer` takes the uppercased reference (`RefSeq::fetch_into`, the canonical
`{A,C,G,T,N}` view) — a deliberate divergence from production, decided 2026-07-26.** An alignment has
no business treating a base differently because a masking tool lowercased it, and ng's alignments work
on uppercase sequences so that no such discrimination is possible.

**What that fixes.** Production hands `F3` the *raw* bytes, and its shift loop compares raw bytes with
no case folding ([norm_seqs.rs](../../../../src/norm_seqs.rs) `next_base_on_left_is_same`), while read
bases are uppercased at decode ([alignment_input.rs](../../../../src/bam/alignment_input.rs)). So on a
soft-masked reference the read shows `A`, the reference shows `a`, they compare unequal, and the shift
loop never advances: **left-alignment silently does nothing in soft-masked regions** — which is to say,
in repeats and low-complexity sequence, exactly where indel placement is ambiguous and left-alignment
is the point. Production's own reason for raw bytes is not a design one — the module note says it
"mirrors exactly what the old serial path did" and is kept for **byte-identity**
([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs), module note) — and the
sibling filter in the same fold, `F1`, *does* fold case
([alignment_input.rs](../../../../src/bam/alignment_input.rs) `read_exceeds_mismatch_fraction`), which
is what marks the left-aligner's case-sensitivity as an oversight rather than an intent.

**Measured, 2026-07-26 — the mechanism is confirmed and the exposure is small.** The parity fixture
(§9) demonstrates it directly: on a soft-masked reference production returns the *mapper's own CIGAR*,
byte for byte, for every read with a shiftable indel, while ng's answer is unchanged from its
unmasked one. But the references actually in use are almost entirely clean, so **this is a latent
robustness fix, not an active defect**:

| reference | lowercase bases | exposure |
|---|---|---|
| GRCh38 no-alt analysis set | **0** (full scan) | none — the GIAB results are unaffected |
| tomato SL4.0 | **227,170** (full scan), across ch01/03/05–12 | ~0.03% of the genome, in callable chromosome sequence |

So it **cannot** explain production's indel deficit, which was measured on GRCh38. (It would bite a
run against UCSC's `hg38.fa` or a typical RepeatMasker-soft-masked plant genome, where the masked
fraction is most of the genome rather than a rounding error.)

One consequence worth naming: canonicalisation also folds non-ACGT to `N`, so an `N` in the reference
now compares *equal* to an `N` in a read where before it did not. A gap may therefore shift across an
`N` run. This is a corner with no practical weight — indels inside `N` runs are junk either way — but
it is a behaviour change, not just a case change.

**ng therefore needs one reference view, where production carries two** (raw for `F3`/`F1`, canonical
for BAQ). Step 1's mismatch filter still fetches raw to stay parity-exact with production's `F1`
([read_filtering.md](read_filtering.md)); whether it should also move to canonical is that spec's
question, not this one. Otherwise this mirrors step 1's `ReadFilter` shape (`reference: R` + a reused
buffer).

**Where it runs — on the worker threads, not on the walker.** Production folded `G2/F3/F1/BAQ` into
one `process_read` and moved it *off* the serial path precisely to shrink the serial floor to the
coordinate merge plus the walker; it is called inside a rayon `map_init` whose per-worker state is the
raw-reference cache and the BAQ engine
([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs) module note;
[baq_stream.rs](../../../../src/pileup/per_sample/baq_stream.rs)). ng follows that shape: preparation
runs per read on the workers, **one `Scratch` per worker**, and the walker only *consumes*
`PreparedRead`s. Read preparation **composes** with the gatherer, it does not fuse into it — which
keeps the bake-off surface alive.

---

## 7. Error model

Two outcomes, never confused:

- **A read produces no usable observation** — `prepare_read` returns `Ok(None)`, the reason tallied, and
  the run continues. **In v1 this never happens on the generic path:** the only step that could decline
  a read was BAQ, which is deferred (§10). The `Option` is kept for the contract and for the modes that
  will reintroduce a decline (BAQ, if it returns; a re-align that cannot place a read).
- **A reference fetch fails** — a contig mismatch, a window past a contig end. A broken run: it is
  **returned as `Err`**, never folded into a per-read `Ok(None)` and never a panic. Returning it is the
  general policy where there is no obvious way to deal with a failure locally, and there is none here;
  a panic would also be the wrong shape on a worker thread, where the driver wants a value it can
  report with the contig and position. This follows step 1, which made the same call and gave its
  reason: a validly-aligned read never covers positions the contig does not have, so an out-of-bounds
  fetch means a malformed record — corrupt input to fail loudly on
  ([read_filtering.md](read_filtering.md) §7). **It is a deliberate divergence from production, and
  therefore from the parity oracle — see §9.** With the fetch now conditional (§5), only indel-bearing
  reads can reach it at all.

Every `Ok(None)` is counted by reason (a per-sample tally, the analogue of `ReadFilterCounts`); v1's
tally is simply all-zero.

---

## 8. Cross-cutting concerns

- **Performance / parallelism.** Pairwise-independence makes preparation embarrassingly parallel with
  deterministic output. v1 is cheap, and cheaper than production's fold: a read with no indel costs one
  CIGAR scan, no reference fetch and no buffer (§5). The re-align mode (gated) would be the per-read
  cost, and it is why the scratch (§6) is reused rather than allocated per read.
- **Determinism.** No mate/column context, so the same read prepares to the same `PreparedRead`
  regardless of thread interleaving — a property downstream determinism relies on.

---

## 9. Reuse over rewrite — the map to production, and the parity oracle

The parity oracle is the production prepared read **with BAQ off**: a ported v1 impl is correct when its
`PreparedRead` is byte-identical to production's `--no-baq` output on a fixture — same canonicalized
CIGAR, and the raw qualities copied through uncapped. Three things to know before writing it:

- **How to run it.** In-process, against `process_read(read, None /* no BAQ */, &mut raw_ref, &cfg)`,
  with **`max_read_mismatch_fraction: None`**. ng moved `F1` to step 1, so leaving it on makes
  production drop reads ng keeps and the two keep-sets diverge for reasons that have nothing to do
  with preparation.
- **What it actually proves — the window fetch and the field wiring, not left-alignment.**
  Left-alignment parity is already banked: ng's default normalizer *is* production's
  `left_align_indels` behind a wrapper, byte-parity asserted where it landed
  ([alignment_normalization.md](../impl_plan/alignment_normalization.md) step B1).
- **Where parity stops, deliberately — two places.**
  - *A failed reference fetch.* Zero span, unknown contig, repository miss, `pos == 0`, a window past
    the contig end: production **skips `F3`/`F1` and emits the read un-left-aligned**
    ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs), `fetch_raw_slice` +
    "Reference-fetch failure is not an error"). ng makes that fatal (§7). The fixture must exclude
    those reads or assert ng's abort; it cannot assert byte-parity on them.
  - *A soft-masked (or non-ACGT) reference.* ng uppercases, production does not (§6), so on lowercase
    reference the two differ **by design** — and this is the useful part: **run the fixture on a
    soft-masked reference deliberately, and the divergence is the measurement.** Every read whose
    indel ng shifts and production leaves put is a read production is mis-placing today. A fixture on
    an all-uppercase reference must show byte-parity; one on a masked reference must not, and the
    difference is the size of the defect.

| what | existing code | status / ng reuse |
|---|---|---|
| indel left-alignment | `AlignmentNormalizer` + `StructuredLeftAligner`/`RepeatedLeftAligner`/`FixpointLeftAligner` ([src/ng/alignment/](../../../../src/ng/alignment/)) | **built** — the preparer calls a normalizer; parity-checked vs production's `left_align_indels` |
| the per-read prep fold | `process_read` ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | model for `LeftAlignPreparer` — its **F3 stage only** (G2/F1 are step-1 filters; BAQ is deferred, §10) |
| building the `PreparedRead` | `prepare_passthrough` (the `--no-baq` path, [baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)) | **`pub`, and it does all of it** — `alignment_end`, `mate_role`, `qname`, `mq_log_err`, adaptor boundary, and the raw-`qual`→`bq_baq` copy. Whether v1 *calls* it or ports it is **open** (§11) |
| ↳ its inner wiring | `mapped_to_prepared` ([baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)) | **private** — not callable from ng; reachable only through `prepare_passthrough` above |
| the prepared read | `PreparedRead` ([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)) | **copy into `src/ng/read/` and extend with `read_group`** — reversed 2026-07-28 (§3 fold-in); production is frozen and the group must survive the preparation boundary. Its misplacement inside `pileup/walker/` is unchanged and still deferred to the port-back |
| the re-align aligner | general/affine best-path aligner ([`alignment.md`](alignment.md) §4.1) | **not built — gated** (algorithm 2, `alignment_best_path.md` Milestone E); its trigger is open too (§4) |
| reference | `RefSeq` — the canonical (uppercased) view ([ref_seq.md](ref_seq.md)) | reuse as-is. **Not** production's raw view: ng uppercases so masking cannot change an alignment (§6). One view, where production carries two |

---

## 10. Deferred / out of scope, with a home

- **BAQ (base alignment quality) — deferred sine die** (owner, 2026-07-25: not of interest now). BAQ is
  a banded HMM (`probaln_glocal`) that estimates, per base, the probability a base is *mis-aligned* and
  caps its quality at that confidence (`bq = min(base_quality, BAQ)`), de-weighting bases in ambiguous
  indels. When/if it returns it is a **config toggle** on the preparer (production's `--no-baq` is the
  off-position, which is v1), not a second implementation — reuse `BaqEngine::process`
  ([baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)). It would reintroduce an
  `Ok(None)` decline reason (§7); it needs no *extra* reference view, since v1 already fetches the
  canonical one (§6). Its home is a later config mode on `LeftAlignPreparer`.
- **The re-align mode (algorithm 2, affine aligner + its trigger) — gated** (§4). Not out of scope, but
  blocked on both a built affine aligner and a not-to-be-trusted judgement.
- **Adaptor-mask application, mate-overlap reconciliation, CIGAR decomposition → the pileup walker.**
  `PreparedRead` is *decomposable*, not decomposed; it *carries* the adaptor boundary but does not apply
  it.
- **Local haplotype reassembly — out of scope, not deferred** (§1).
- **freebayes' iterated `stablyLeftAlign` → a refinement.** freebayes iterates `leftAlign` to
  convergence; `FixpointLeftAligner` is ng's counterpart. Whether that is the whole of freebayes'
  indel-stabilisation needs pinning against `freebayes/src/LeftAlign.cpp` before implementing, then
  measuring whether it moves indel placement on our data. (freebayes is MIT — reading and porting is
  fine.)

---

## 11. Resolved decisions & open questions

**Resolved.**

- **Read preparation is a generic-path-only step (2026-07-25).** The STR path re-aligns every spanning
  read against the tract, so canonicalizing the mapper's CIGAR there would be work nothing reads — and
  would perturb the slice the re-alignment is given (§1). Its per-read operation also produces an
  observation, not a read, so it belongs to observation generation
  ([`locus_generation_ssr.md`](locus_generation_ssr.md)). The two path specs are retired and the
  `ReadPreparer` trait loses its path-owned associated types (§6).
- **BAQ deferred sine die (2026-07-25).** So v1 is pass-through + left-alignment, `bq_baq` is raw
  qualities uncapped, the parity oracle is production's `--no-baq`, and the generic preparer never
  returns `Ok(None)` (§7). The toggle design is recorded for if BAQ returns (§10).
- **The re-align mode is gated** — its affine aligner (algorithm 2) is unbuilt and its trigger is open
  (§4). Not out of scope; simply not v1.
- **~~Reuse production's `PreparedRead`~~ — REVERSED 2026-07-28: ng owns its own copy, extended with `read_group`** (§3 fold-in). The rest stands: the reference is held, not passed; dispatch is static. (§3, §6.)
- **The signature (2026-07-25):** `prepare_read(&self, read: MappedRead, scratch) -> Result<Option<PreparedRead>, ReadPrepError>`.
  By value; `Result` around `Option` because a decline and a broken run are different outcomes; the
  reference held per impl rather than an `Option<R>` on one type. Rationale and the alternatives beaten
  are in §6.
- **Left-alignment is the reference's only consumer, so the fetch is conditional on the read carrying
  an indel (2026-07-25).** Production fetches for every read because `F1` needed the same window; ng
  moved `F1` to step 1, so that reason does not carry over (§5).
- **Alignments work on uppercase sequences (owner, 2026-07-26).** There is no reason for an alignment
  to treat a base differently because a masking tool lowercased it, so `LeftAlignPreparer` fetches the
  canonical (uppercased) reference rather than production's raw bytes (§6). Production's raw choice was
  never designed — it is kept for byte-identity with an older serial path — and it has a consequence:
  because read bases are uppercased at decode while the reference is not, production's left-alignment
  **does nothing on a soft-masked reference**. That is recorded here as a **known production defect**
  for the port-back; production is frozen, so ng fixes it on its own side and §9's fixture measures it.
  **It is latent, not active** (§6): GRCh38 carries no lowercase at all and tomato only 227 kb, so it
  explains nothing we have observed — including the indel deficit, which was measured on GRCh38.
- **Reusing production: call it where it fits, decided case by case (owner, 2026-07-25).** ng may call
  a production function when it fits ng's need as-is — the default normalizer already does
  ([left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs) wraps
  `left_align_indels`). What is **not** allowed is contorting ng's code to avoid a port. So each reuse
  is judged on its own; there is no blanket rule either way.

**Open (confirm before the relevant code).**

- **What flags a region as "not to be trusted", and when?** (§4) The re-align trigger — a new producer
  or a two-pass arrangement. Settle (together with building algorithm 2) before that mode ships; it
  decides whether the mode is worth its cost. **This step now carries the affine aligner's cost:** the
  alignment plan gates Milestone E on this question and hands the milestone over if the gate is still
  open at its Checkpoint D ([alignment_best_path.md](../impl_plan/alignment_best_path.md) Milestone E).
  [`../impl_plan/read_preparation.md`](../impl_plan/read_preparation.md) records it as a follow-on
  plan, since v1 does not build the mode.
- **Call production's `prepare_passthrough`, or port it?** (§9) *Leaning: call it in the first slice;
  port if the copy below shows up in a profile.* The policy is settled (see Resolved above) and this
  call fits it — both types are already ng's, so nothing bends. The one cost found by reading it:
  `prepare_passthrough` does `read.qual.clone()` while `mapped_to_prepared` never reads `read.qual`, so
  calling it pays an avoidable `Vec<u8>` copy per read (production pays it for symmetry with the BAQ
  arm, where `bq_baq` genuinely differs); a port would `mem::take` instead. Against that, calling makes
  the field wiring parity-exact by construction. Confirm before the code.

*(The former "does pass-through skip the base-quality capping too?" question is moot — there is no
capping step in v1. The former STR-path open questions moved with the STR path to
[`locus_generation_ssr.md`](locus_generation_ssr.md).)*
