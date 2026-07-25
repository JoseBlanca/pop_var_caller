# ng — read preparation (the generic SNP/indel path)

*Status: design spec, rewritten 2026-07-25. **This is the single read-preparation spec** — read
preparation is a **generic-path-only** step (§1 says why), so there is no shared preamble and no STR
sibling; the former `read_preparation_generic.md` and `read_preparation_ssr.md` are short redirects
here. **What is built, and what v1 is:** the left-alignment transform ships (`AlignmentNormalizer` +
three impls in `src/ng/alignment/`); **BAQ is deferred sine die** (§10 — not of interest now); the
**re-align mode's aligner (algorithm 2, affine) and its trigger are both unbuilt/gated** (§4, §10); and
the **`PreparedRead`-producing preparer is not built**. So **v1 read preparation is pass-through +
left-alignment** — nothing else. Grounded in the production `process_read` fold
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
not generate candidate alleles or compute read likelihoods (step 7). It does not recalibrate base
qualities (BQSR). And **it does not reassemble: local haplotype reassembly is out of scope for ng, not
deferred** — the production caller already calls generic loci better than GATK without it, and it would
break the per-read independence every mode below relies on (§6).

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
indels has nothing to shift, so canonicalizing it is provably a no-op. Recognizing that from the read's
own alignment record and skipping the work changes nothing. (With BAQ deferred, §10, pass-through is the
whole of "do nothing" — there is no quality-capping step left to decide about.)

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
transform's output, field for field, exactly as step 1 reuses `MappedRead`. Reproduced for the reader,
**not redefined**:

```rust
pub struct PreparedRead {
    pub chrom_id: ContigId,
    pub alignment_start: Position,
    pub alignment_end: Position,          // cached, so the walker never re-walks the CIGAR for span
    pub cigar: Vec<CigarOp>,              // canonicalized (or, later, re-aligned) — unlike MappedRead.cigar
    pub seq: Vec<u8>,                     // uppercase ACGTN
    pub bq_baq: Vec<u8>,                  // in v1: the read's base qualities, passed through UNCAPPED
                                          //   (BAQ deferred, §10). The field is production's; a future
                                          //   BAQ would cap it to min(BQ, BAQ).
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
not. That needs nothing but the read, and it is all of v1.

**Re-align is different: it needs a judgement the read cannot make about itself.** "The mapper's answer
here is not trustworthy" is a property of the *place*, not one read — a region where reads disagree,
pile up mismatches in one column, or clip at the same offset. Nothing in the current ng step map
produces that judgement:

- **Region typing** classifies the reference (microsatellite, cluster, satellite, generic). It says
  what the reference *is*, not how well reads mapped to it — it never looks at reads.
- **The evidence-gatherer** does see the reads and could discover it, but it runs *after* preparation,
  so a verdict it produces arrives too late for the read it should have changed.

So the trigger needs either a new producer or a deliberate two-pass arrangement, and picking one is open
(§11). It matters because **the trigger, not the algorithm, sets how much this mode is worth**: an
aligner that never fires rescues nothing. This is one of the two reasons the re-align mode is gated; the
other is that its aligner (algorithm 2) is unbuilt.

---

## 5. The transform detail

Production runs one per-read fold, `process_read`, whose stages are
`G2 bad-CIGAR → F3 left-align → F1 mismatch-fraction → BAQ`. ng assigns the two *rejects* — `G2` and
`F1` — to **step 1** (filters #9 and #8 in `read_filtering.md` §3), and **defers `BAQ` sine die**
(§10). What remains for v1 preparation is `F3`:

- **Indel left-alignment** (the canonicalize mode) — call an `AlignmentNormalizer` (§2). Rewrites
  **only the CIGAR**; bases and qualities are untouched. **Built.** `bq_baq` is filled by copying the
  read's raw qualities through uncapped (§3), production's `--no-baq` behaviour.

**Why the step-1/step-2 split is safe.** Production runs `F1` *after* left-alignment; ng runs mismatch
filtering in step 1, *before* preparation's left-alignment. Safe because **left-alignment provably
preserves the mismatch count** — a debug-assert in production's `left_align_indels` guarantees it — so
ng's order gives the identical verdict, and the bad-CIGAR check (`G2`/#9) sees the raw decoded CIGAR in
both.

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
read-likelihood models do — it holds the normalizer's working buffers now, and the alignment matrices
the re-align mode will need later.

```rust
pub trait ReadPreparer {
    /// Reused buffers — allocated once per worker, never per read.
    type Scratch: Default;
    /// Prepare one filtered read against the reference around its OWN span. No locus argument
    /// (the generic transform is locus-independent) and no window argument — the impl HOLDS its
    /// reference accessor and fetches what it needs. `None` if unusable (tallied); see §7 — in v1
    /// nothing returns `None`.
    fn prepare_read(&self, read: &MappedRead, scratch: &mut Self::Scratch) -> Option<PreparedRead>;
}
```

The v1 implementation is **`LeftAlignPreparer`** (pass-through + canonicalize). The trait carries **no
`type Locus`** (the generic path needs none) and **no `type Prepared`** (always `PreparedRead`); the
path-owned associated types an earlier sketch had — which could not compile as a `Box<dyn ReadPreparer>`
— are gone with the STR arm. The trait still exists because a second impl (the re-align mode) will bench
behind it.

The implementation **holds its own reference accessor** as a field and fetches around each read's span;
there is no reference-window argument. v1 needs only **raw, case-preserving bytes** (`RawRefSeq`) for
left-alignment (the aligner's own view). (A **canonical** uppercased view was needed only by BAQ, which
is deferred — §10 — so v1 touches only the raw accessor.) This mirrors step 1's `ReadFilter`
(`reference: R` + a reused `ref_buf`).

**Where it is invoked — by the pileup, per read, as it walks a non-STR stretch.** Following production,
where `process_read` runs as the walker ingests reads: the pileup **calls** `prepare_read` on each read
and consumes the `PreparedRead`. Read preparation **composes** with the gatherer, it does not fuse into
it — which keeps the bake-off surface alive.

---

## 7. Error model

Two outcomes, never confused:

- **A read produces no usable observation** — `prepare_read` returns `None`, the reason tallied. **In v1
  this never happens on the generic path:** the only step that could decline a read was BAQ, which is
  deferred (§10). The `Option` is kept for the contract and for the modes that will reintroduce a
  decline (BAQ, if it returns; a re-align that cannot place a read).
- **A reference fetch fails** — a contig mismatch, a window past a contig end. A broken run, **fatal**,
  surfaced as such, never folded into a per-read `None`.

Every `None` is counted by reason (a per-sample tally, the analogue of `ReadFilterCounts`); v1's tally is
simply all-zero.

---

## 8. Cross-cutting concerns

- **Performance / parallelism.** Pairwise-independence makes preparation embarrassingly parallel with
  deterministic output. v1 (left-alignment) is cheap; the re-align mode (gated) would be the per-read
  cost, and it is why the scratch (§6) is reused rather than allocated per read.
- **Determinism.** No mate/column context, so the same read prepares to the same `PreparedRead`
  regardless of thread interleaving — a property downstream determinism relies on.

---

## 9. Reuse over rewrite — the map to production, and the parity oracle

The parity oracle is the production prepared read **with BAQ off**: a ported v1 impl is correct when its
`PreparedRead` is byte-identical to production's `--no-baq` output on a fixture — same canonicalized
CIGAR, and the raw qualities copied through uncapped. One parity fixture, which also proves left-alignment
in isolation, since nothing else touches the read.

| what | existing code | status / ng reuse |
|---|---|---|
| indel left-alignment | `AlignmentNormalizer` + `StructuredLeftAligner`/`RepeatedLeftAligner`/`FixpointLeftAligner` ([src/ng/alignment/](../../../../src/ng/alignment/)) | **built** — the preparer calls a normalizer; parity-checked vs production's `left_align_indels` |
| the per-read prep fold | `process_read` ([read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | model for `LeftAlignPreparer` — its **F3 stage only** (G2/F1 are step-1 filters; BAQ is deferred, §10) |
| `bq_baq` passthrough | `prepare_passthrough` (the `--no-baq` path, [read_processor.rs](../../../../src/pileup/per_sample/read_processor.rs)) | model — copies raw `qual` into `bq_baq` uncapped; this **is** v1's quality handling |
| the prepared read | `PreparedRead` + `mapped_to_prepared` ([pileup/walker/mod.rs](../../../../src/pileup/walker/mod.rs)) | **reuse as-is**; may want hoisting out of `pileup/walker/` |
| the re-align aligner | general/affine best-path aligner ([`alignment.md`](alignment.md) §4.1) | **not built — gated** (algorithm 2, `alignment_best_path.md` Milestone E); its trigger is open too (§4) |
| reference | `RawRefSeq` ([ref_seq.md](ref_seq.md)) | reuse as-is (raw only in v1; canonical was BAQ's, deferred) |

---

## 10. Deferred / out of scope, with a home

- **BAQ (base alignment quality) — deferred sine die** (owner, 2026-07-25: not of interest now). BAQ is
  a banded HMM (`probaln_glocal`) that estimates, per base, the probability a base is *mis-aligned* and
  caps its quality at that confidence (`bq = min(base_quality, BAQ)`), de-weighting bases in ambiguous
  indels. When/if it returns it is a **config toggle** on the preparer (production's `--no-baq` is the
  off-position, which is v1), not a second implementation — reuse `BaqEngine::process`
  ([baq_engine.rs](../../../../src/pileup/per_sample/baq_engine.rs)). It would reintroduce a `None`
  decline reason (§7) and require the **canonical** reference view (§6). Its home is a later config mode
  on `LeftAlignPreparer`.
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

- **Read preparation is a generic-path-only step (2026-07-25).** The locus-independence test (§1): the
  STR "prep" needs the locus and produces an observation, not a read, so it is observation generation
  ([`locus_generation_ssr.md`](locus_generation_ssr.md)). The two path specs are retired and the
  `ReadPreparer` trait loses its path-owned associated types (§6).
- **BAQ deferred sine die (2026-07-25).** So v1 is pass-through + left-alignment, `bq_baq` is raw
  qualities uncapped, the parity oracle is production's `--no-baq`, and the generic preparer never
  returns `None` (§7). The toggle design is recorded for if BAQ returns (§10).
- **The re-align mode is gated** — its affine aligner (algorithm 2) is unbuilt and its trigger is open
  (§4). Not out of scope; simply not v1.
- **Reuse production's `PreparedRead`; the reference is held, not passed; dispatch is static.** (§3, §6.)

**Open (confirm before the relevant code).**

- **What flags a region as "not to be trusted", and when?** (§4) The re-align trigger — a new producer
  or a two-pass arrangement. Settle (together with building algorithm 2) before that mode ships; it
  decides whether the mode is worth its cost.

*(The former "does pass-through skip the base-quality capping too?" question is moot — there is no
capping step in v1. The former STR-path open questions moved with the STR path to
[`locus_generation_ssr.md`](locus_generation_ssr.md).)*
