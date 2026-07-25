# ng — read preparation (the generic SNP/indel path)

*Status: design spec, rewritten 2026-07-25. **This is the single read-preparation spec** — read
preparation is a **generic-path-only** step (§1 says why), so there is no shared preamble and no STR
sibling any more; the former `read_preparation_generic.md` and `read_preparation_ssr.md` are now short
redirects here. Defines what read preparation does, how much work it gives each read, and how it
composes with the alignment module ([`alignment.md`](alignment.md)). **Partly built:** the
left-alignment transform ships (`AlignmentNormalizer` + three impls in `src/ng/alignment/`) and so does
the re-align mode's best-path aligner; **BAQ, the re-align *trigger*, and the `PreparedRead`-producing
preparer are not built yet** (§9). Grounded in the production `process_read` fold
([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)). Naming: **STR** in prose,
`ssr` in code.*

---

## 1. What read preparation is — and why only the generic path has it

Read filtering decided *which* reads survive ([`read_filtering.md`](read_filtering.md)). Read
preparation is a **per-read, locus-independent transform**: it canonicalises a filtered `MappedRead`
*against the reference around the read's own span*, producing a `PreparedRead` that **is still a read** —
decomposable, and reused across *every* locus it overlaps. It runs before, and independent of, which
loci exist.

**That locus-independence is the definition, and it is what makes read preparation a generic-path-only
step.** The two marker paths differ exactly here:

- **Generic (SNP/indel).** The transform works against the reference around the read's own span — no
  locus needed. The result serves every position the read covers. *That* is preparing a read.
- **STR.** The per-read operation aligns the read *against a specific tract* to read out what that read
  shows *about that locus*. It needs the locus (motif, borders, flanks), and its output is not a
  canonicalised read but an **observation** — the same read at another tract comes out differently. That
  is not preparing a read; it is **generating a locus's evidence from a read**. It lives in the STR
  observation generator ([`locus_generation_ssr.md`](locus_generation_ssr.md)), which calls the
  repeat-aware aligner ([`alignment.md`](alignment.md) §4.2) per read. So **the STR path has no read
  preparation** — it goes filtering → observation generation, with the alignment as a component of the
  latter.

The test is simply: *does the transform need the locus to interpret the read?* No → preparation. Yes →
observation generation. An earlier design forced both paths under one `ReadPreparer` trait with a
path-owned output and a "the STR path always aligns" mode; that unified two operations that are not the
same kind, and it is retired.

**Non-goals.** Preparation never *drops* a read for a whole-read property — that was filtering. It does
not decompose a read into per-position events, apply the adaptor mask, or reconcile overlapping-mate
qualities — those need the locus-column context and are the **pileup walker's** job (§8, §10). It does
not generate candidate alleles or compute read likelihoods (step 7). And **it does not reassemble:
local haplotype reassembly is out of scope for ng, not deferred** — the production caller already calls
generic loci better than GATK without it, so it buys nothing here, and it would break the per-read
independence every mode below relies on (§7).

---

## 2. How much work a read needs — the three modes

Preparation picks one of three modes per read. All three produce the **same output type**
(`PreparedRead`, §3), which is what makes them interchangeable and comparable.

| mode | what it does | when | cost |
|---|---|---|---|
| **pass through** | nothing to the placement | the read shows no insertions or deletions | none |
| **canonicalize** | rewrite the read's indels into their leftmost equivalent spelling; optionally cap base qualities by alignment confidence | the read has indels and its placement is trusted | cheap |
| **re-align** | discard the mapper's line-up and compute a fresh one from the read's bases | the read's placement is not trusted (§4) | a full alignment per read |

**Pass-through is a fast path, not a different answer.** Left-alignment shifts indels; a read with no
indels has nothing to shift, so canonicalizing it is provably a no-op. Recognizing that from the read's
own alignment record and skipping the work changes nothing. (Whether it also skips the base-quality
capping — which reads alignment *confidence*, not indel placement — is open, §9.)

**Canonicalize is about spelling, not quality.** The same indel can be written at several equivalent
reference positions when it sits in or near a repeat — the gap slides without changing a base of the
result. Left-alignment picks the leftmost spelling so equivalent variants get an identical one;
otherwise the reads supporting them scatter across several weak candidates instead of pooling into one
strong one. The operation lives in the alignment module ([`alignment.md`](alignment.md) §6) as the
`AlignmentNormalizer` trait (**built** — `StructuredLeftAligner`, the GATK/production port and default;
`RepeatedLeftAligner`, the freebayes-shaped repeated-pass form; `FixpointLeftAligner`, the fail-loud
fixpoint wrapper). Preparation *calls* a normalizer; it does not re-implement left-alignment. The
optional quality cap is **BAQ** (§5).

**Re-align is the only mode that questions the mapper.** The other two accept the read's placement; this
one throws it away and computes a new line-up with a best-path alignment algorithm
([`alignment.md`](alignment.md) §4.1). It is the expensive, rare mode, and the only route by which a
mis-placed read is rescued — so its trigger (§4) sets how much of that class the caller ever recovers.

---

## 3. The output — `PreparedRead` (reused from production)

Every mode yields the **existing production `PreparedRead`**
([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)), **reused as-is** — it already *is* the
transform's output, field for field, exactly as step 1 reuses `MappedRead`. Reproduced for the reader,
**not redefined**:

```rust
pub struct PreparedRead {
    pub chrom_id: ContigId,
    pub alignment_start: Position,
    pub alignment_end: Position,          // cached, so the walker never re-walks the CIGAR for span
    pub cigar: Vec<CigarOp>,              // canonicalized or re-aligned (unlike MappedRead.cigar)
    pub seq: Vec<u8>,                     // uppercase ACGTN
    pub bq_baq: Vec<u8>,                  // BAQ-capped min(BQ, BAQ) (unlike MappedRead.qual)
    pub mq_log_err: f64,                  // derived ln(P_err) from MAPQ (not on MappedRead)
    pub mapq: MapQual,                    // raw, preserved
    pub is_reverse_strand: bool,          // decoded from flag
    pub mate_role: MateRole,              // Solo | FirstOfPair | SecondOfPair (from flag bits)
    pub adaptor_boundary: Option<Position>, // carried through, applied later by the walker (§10)
}
```

A `PreparedRead` carries **no** per-base overlap adjustment — that is pairwise and happens in the walker
(§10), which is what keeps it self-contained and pairwise-independent. Reused unchanged; the walker port
may want it hoisted out of `pileup/walker/`, since preparation produces it and the pileup only consumes
it.

---

## 4. Choosing the mode — the part that is not settled

Pass-through and canonicalize are chosen from the read's own alignment record: does it carry indels or
not. That needs nothing but the read.

**Re-align is different: it needs a judgement the read cannot make about itself.** "The mapper's answer
here is not trustworthy" is a property of the *place*, not one read — a region where reads disagree,
pile up mismatches in one column, or clip at the same offset. Nothing in the current ng step map
produces that judgement:

- **Region typing** classifies the reference (microsatellite, cluster, satellite, generic). It says
  what the reference *is*, not how well reads mapped to it — it never looks at reads.
- **The evidence-gatherer** does see the reads and could discover it, but it runs *after* preparation,
  so a verdict it produces arrives too late for the read it should have changed.

So the trigger needs either a new producer or a deliberate two-pass arrangement, and picking one is open
(§9). It matters because **the trigger, not the algorithm, sets how much this mode is worth**: an
aligner that never fires rescues nothing.

---

## 5. The transform detail — left-align, BAQ, re-align

Production runs one per-read fold, `process_read`, whose stages are
`G2 bad-CIGAR → F3 left-align → F1 mismatch-fraction → BAQ`. ng already assigned the two *rejects* —
`G2` and `F1` — to **step 1** (filters #9 and #8 in `read_filtering.md` §3). What remains for
preparation is the transforms:

1. **Indel left-alignment** (the canonicalize mode) — as §2: call an `AlignmentNormalizer`. Rewrites
   **only the CIGAR**; bases and qualities are untouched. **Built.**
2. **BAQ (base alignment quality) — optional.** A banded HMM re-aligns the read to the reference to
   estimate, per base, the probability it is *mis-aligned*, and caps each base quality at that
   confidence (`bq = min(base_quality, BAQ)`); bases the HMM cannot place are set to quality 0. BAQ
   de-weights bases in and near ambiguous indels **without rewriting the alignment**. It can decline a
   read outright (HMM overflow; reference window past the contig end; no aligned `M` op) → `None` (§7).
   An htslib `probaln_glocal` port. **Not built yet.**
3. **Re-alignment** (the re-align mode) — a best-path alignment from `alignment.md` §4.1 that replaces
   the CIGAR wholesale. **The aligner is built; its trigger is not** (§4).

**Jargon, once — BAQ.** *Base Alignment Quality* is a per-base confidence that a base is aligned to the
right reference position (as opposed to the base-*call* quality, confidence in the letter). A base in an
ambiguously-placed indel gets a low BAQ and is de-weighted, so a mis-alignment cannot masquerade as a
confident mismatch.

**BAQ is a toggle, not a second implementation.** Left-align-only (BAQ off) is a first-class v1 mode —
the freebayes-style preparation — expressed as a config toggle (`baq: Option<BaqConfig>`), not a sibling
type. Production models it this way: a `--no-baq` path whose `prepare_passthrough` copies raw qualities
into `bq_baq` uncapped. *Rejected alternative:* two sibling impls — they differ *only* by whether the
final cap is applied, so a second preparer buys one bake-off row at the price of two code paths for one
algorithm. Same convention as step 1's `max_read_mismatch_fraction: None`.

**Why the step-1/step-2 split is safe.** Production runs `F1` *after* left-alignment; ng runs mismatch
filtering in step 1, *before* preparation's left-alignment. Safe because **left-alignment provably
preserves the mismatch count** — a debug-assert in production's `left_align_indels` guarantees it — so
ng's order gives the identical verdict, and the bad-CIGAR check (`G2`/#9) sees the raw decoded CIGAR in
both.

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

**Buffers are caller-owned and reused.** The alignment algorithms preparation calls need matrices;
allocating them per read is the other cost this path cannot carry. So preparation threads a reusable
scratch value, as the read-likelihood models do.

```rust
pub trait ReadPreparer {
    /// Reused buffers, including those the alignment algorithms need — allocated once per
    /// worker, never per read.
    type Scratch: Default;

    /// Prepare one filtered read against the reference around its own span. `None` = no usable
    /// observation here (§7), tallied. The impl HOLDS its own reference accessors as fields.
    fn prepare_read(&self, read: &MappedRead, scratch: &mut Self::Scratch) -> Option<PreparedRead>;
}
```

The trait exists because there is a **second generic implementation to bench behind it** — the re-align
mode's is one axis, and a swap of the best-path algorithm inside `alignment.md` is another — not to
abstract over paths. So it carries **no `type Locus`** (the generic path needs none) and **no
`type Prepared`** (always `PreparedRead`); the path-owned associated types an earlier sketch carried —
which could not even compile as a `Box<dyn ReadPreparer>` — are gone with the STR arm.

The implementation **holds its own reference accessors** as fields and fetches around each read's span;
there is no reference-window argument. The transform needs two views of the reference at once — **raw,
case-preserving bytes** (`RawRefSeq`) for left-alignment (the aligner's own view) and **canonical
uppercased bytes** (`RefSeq`) for the BAQ HMM — which a single materialised window cannot carry. With
BAQ off, only the raw accessor is touched. This mirrors step 1's `ReadFilter` (`reference: R` + a reused
`ref_buf`).

**Where it is invoked — by the pileup, per read, as it walks a non-STR stretch.** Following production,
where `process_read` runs as the walker ingests reads: the pileup **calls** `prepare_read` on each read
and consumes the `PreparedRead`. Read preparation **composes** with the gatherer, it does not fuse into
it — which keeps the bake-off surface alive (a re-aligning preparer, or a different alignment algorithm
behind one, swaps in behind the same contract).

---

## 7. Error model

Two outcomes, never confused:

- **A read produces no usable observation** — `prepare_read` returns `None`, the reason tallied. Normal.
  On the generic path the reason is `Baq` (BAQ declined the read — HMM overflow, reference window past
  the contig end, no aligned `M` op, an `N`/ref-skip in the CIGAR). With BAQ off, the transform never
  returns `None`.
- **A reference fetch fails** — a contig mismatch, a window past a contig end. A broken run, **fatal**,
  surfaced as such, never folded into a per-read `None`. So a `None` always means "this read, unusable
  here" and never hides a broken reference.

Every `None` is counted by reason (a per-sample tally, the analogue of `ReadFilterCounts`). A run that
silently prepares nothing must be distinguishable from one that prepared everything only by the counts.

---

## 8. Cross-cutting concerns

- **Performance / parallelism.** Pairwise-independence makes preparation embarrassingly parallel with
  deterministic output. BAQ (the per-read HMM) is the cost; reuse the per-worker scratch (§6) rather
  than allocating per read.
- **Determinism.** No mate/column context, so the same read prepares to the same `PreparedRead`
  regardless of thread interleaving — a property downstream determinism relies on.

---

## 9. Reuse over rewrite — the map to production, and the parity oracle

The parity oracle is the production prepared read, **in both quality modes**: a ported impl is correct
when its `PreparedRead` is byte-identical to production's on a fixture — same canonicalized CIGAR, and
same qualities under **BAQ on** (vs production's default) *and* **BAQ off** (vs `--no-baq`). Two parity
fixtures; the BAQ-off one also proves left-alignment in isolation, since nothing else touches the
qualities.

| what | existing code | status / ng reuse |
|---|---|---|
| indel left-alignment | `AlignmentNormalizer` + `StructuredLeftAligner`/`RepeatedLeftAligner`/`FixpointLeftAligner` ([src/ng/alignment/](../../../../src/ng/alignment/)) | **built** — the preparer calls a normalizer; parity-checked vs production's `left_align_indels` |
| the re-align aligner | best-path aligner ([`alignment.md`](alignment.md) §4.1, `src/ng/alignment/`) | **built** — call it; only the *trigger* (§4) is open |
| the per-read prep fold | `process_read` ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | model for the preparer — its F3 + BAQ stages only (G2/F1 are step-1 filters) |
| BAQ (on) | `BaqEngine::process` ([baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)), htslib `probaln_glocal` | **to build** — call directly; `None` on BAQ-skip |
| BAQ (off) — the left-align-only mode | `prepare_passthrough` (the `--no-baq` path, same file) | **to build** — copies raw `qual` into `bq_baq` uncapped |
| the prepared read | `PreparedRead` + `mapped_to_prepared` ([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)) | **reuse as-is**; may want hoisting out of `pileup/walker/` |
| reference | `RawRefSeq` (left-align) + `RefSeq` (BAQ) ([ref_seq.md](ref_seq.md)) | reuse as-is |

---

## 10. Deferred / out of scope, with a home

- **Adaptor-mask application, mate-overlap reconciliation, CIGAR decomposition → the pileup walker.**
  `PreparedRead` is *decomposable*, not decomposed; it *carries* the adaptor boundary but does not apply
  it.
- **Local haplotype reassembly — out of scope, not deferred** (§1): production beats GATK on generic
  loci without it, and it breaks per-read independence.
- **freebayes' iterated `stablyLeftAlign` → a refinement.** freebayes iterates `leftAlign` to
  convergence; `FixpointLeftAligner` is ng's counterpart. Whether that is the whole of freebayes'
  indel-stabilisation needs pinning against `freebayes/src/LeftAlign.cpp` before implementing, then
  measuring whether it moves indel placement on our data. (freebayes is MIT — reading and porting is
  fine.)
- **A faster best-path core** (wavefront / difference-recurrence SIMD) for the re-align mode — a swap
  *inside* `alignment.md` §4, invisible here, since every mode yields the same `PreparedRead`.

---

## 11. Resolved decisions & open questions

**Resolved.**

- **Read preparation is a generic-path-only step (2026-07-25).** The locus-independence test (§1) shows
  the STR "prep" needs the locus and produces an observation, not a read: it is observation generation
  ([`locus_generation_ssr.md`](locus_generation_ssr.md)), not preparation. So the two path specs are
  retired and the `ReadPreparer` trait loses its path-owned associated types (§6). *Beaten:* the
  one-trait-two-paths unification (a "the STR path always aligns" mode), which abstracted over two
  operations that are not the same kind and could not compile as `Box<dyn ReadPreparer>`.
- **Reuse production's `PreparedRead`; the reference is held, not passed; dispatch is static; BAQ is a
  toggle; the step-1/step-2 split is safe.** (§3, §5, §6.)

**Open (confirm before the relevant code).**

- **What flags a region as "not to be trusted", and when?** (§4) The re-align trigger — a new producer
  or a two-pass arrangement. Settle before that mode is built; it decides whether the mode is worth its
  cost.
- **Does pass-through skip the base-quality capping too, or only the left-alignment?** (§2) Cheap to
  settle by comparing output on indel-free reads with the capping on and off.

*(The former STR-path open questions — a fast path for clean microsatellite reads, and which reads reach
the STR aligner / the lower-bound selection — moved with the STR path to
[`locus_generation_ssr.md`](locus_generation_ssr.md), where they belong.)*
