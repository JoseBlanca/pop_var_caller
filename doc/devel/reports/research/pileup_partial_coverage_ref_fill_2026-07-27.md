# Production's pileup credits reads with reference bases they never sequenced

*Research note, 2026-07-27. Found while specifying ng's generic locus generator
([../../ng/spec/locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md) §6), which is
a port of the module described here. **Production is frozen: this is a record, not a change.**
Nothing here is measured yet — the size of the effect is an open number, and §6 of that spec makes
counting it an acceptance criterion for the port.*

---

## 1. The behaviour

An open pileup record can span several reference bases — a 5-base deletion creates a 6-base
footprint. In the fold ([open_record.rs:717-824](../../../../src/pileup/walker/open_record.rs#L717))
**every** contributor at the current walker position is folded into **every** affected record, and
the read's haplotype string over the record's footprint is built by `apply_events_to_ref_into`
([:480-574](../../../../src/pileup/walker/open_record.rs#L480)), which fills any part of the
footprint the read's own events do not cover with bases taken from the **reference** —
[:522-531](../../../../src/pileup/walker/open_record.rs#L522) for the head and interior gaps,
[:568-573](../../../../src/pileup/walker/open_record.rs#L568) for the tail.

So a read whose alignment covers only the first two bases of a six-base record contributes a
six-base haplotype = *its two observed bases* + *four reference bases it never sequenced*, and is
then credited as full support (`num_obs: 1`, plus `q_sum`, `fwd`, `mapq_sum`, `mapq_sum_sq`) for
whichever bucket that string lands in — normally REF.

**There is no guard, at any layer.** The only test in the fold is
[:775-778](../../../../src/pileup/walker/open_record.rs#L775):

```rust
if window_events.is_empty() {
    continue;
}
```

whose comment reads *"No events overlapping = the read doesn't observe this record's REF stretch
**at all** — it shouldn't fold."* — one overlapping base admits the read. Record selection is
`find_overlapping` ([:303](../../../../src/pileup/walker/open_record.rs#L303)), documented as
*"non-empty interval intersection"*. And `apply_events_to_ref_into` could not check even if it
wanted to: its parameters are `(allele_seq, record_pos, ref_seq, events)` — **the read's
`alignment_start`/`alignment_end` are not in scope.** The fold does hold `contrib.alignment_start`,
and uses it three lines later for `placed_left`/`placed_start`
([:790-791](../../../../src/pileup/walker/open_record.rs#L790)) — bias counters that nothing in
`src/` reads — but never for a span test.

**The contribution survives the read.** `expire_passed`
([active_read_set.rs:149-176](../../../../src/pileup/walker/active_read_set.rs#L149)) removes the
read from the active set and touches no open record; the record stays open until its footprint end
passes; `finalise` iterates `folded_reads` unconditionally. Once summed, the contribution cannot be
withdrawn — `subtract_contribution` is reachable only from the re-fold path, which needs the read to
still be live.

## 2. `widen` makes it retroactive

When a record grows, `widen` ([:348-419](../../../../src/pileup/walker/open_record.rs#L348)) appends
the new reference bases to **every** allele bucket — including buckets holding reads that have
already left the active set and can never be re-folded. Such a read's observation is extended with
reference bases the walker had not even fetched when the read was folded. `max_record_span` defaults
to 5000, so a 150 bp read can end up credited with support for a haplotype thousands of bases long.

`widen`'s own comment says the append "reproduces exactly what `apply_events_to_ref` would emit if
we re-folded the read against the wider `ref_seq`" — which is true, and is the problem: reproducing
the fill faithfully reproduces the fabrication.

## 3. Why does it happen? — nobody decided it, and the spec says the opposite

`doc/devel/specs/pileup_walker.md:1073-1092`, §"Step 6: REF bucket maintenance":

> After folding all events, the walker checks every `ActiveRead` that **overlaps** `A` but
> contributed no event in this open record's span (i.e. the read's haplotype string under `A` equals
> `A.ref_seq`). Each such read contributes one observation to the REF bucket.
>
> This is what makes `num_obs` for REF correct: REF is "reads that match the reference across the
> open record's **full span**", not "reads that just happen to match at the anchor position."

The property asserted in the second paragraph is the one you would want. The admission rule
prescribed in the first is `overlaps`. And the parenthetical treats *"produced the reference
string"* and *"matched the reference over the full span"* as the same fact — which is exactly the
conflation, since the reference string is manufactured for the positions the read is missing. The
same full-span framing recurs at [:156-163](../../specs/pileup_walker.md).

**Searched and found nothing.** Case-insensitive search of the 1,442-line
`doc/devel/specs/pileup_walker.md` for `partial` / `partially` / `does not span` / `overhang` /
`censor` / `lower bound` / `truncat` / `fabricat` / `witness`: **zero hits.** Same terms across the
seven pileup review reports (`pileup_2026-05-06/09/11`, the freebayes/GATK/samtools comparisons, the
clone audit): the only hits on `apply_events_to_ref` concern tied-anchor sort order (Mi9) and a
memcpy collapse (L11). **No review ever questioned the reference fill.**

So the answer to "why" is: the shape *one haplotype string per read per record* does this when
nothing checks the read's extent, and the spec recorded the intended invariant as though it had been
implemented.

## 4. What the three reference callers do — production is alone

All three are vendored read-only in the sibling checkout `pop_var_caller/`.

**freebayes (MIT) — an explicit span gate, plus a down-weighted partial channel.** This matters most
because production's pileup is freebayes-derived. `getCompleteObservationsOfHaplotype`
(`freebayes/src/AlleleParser.cpp:3584`):

```cpp
if (ra.start <= currentPosition && ra.end >= currentPosition + haplotypeLength) {
```

Start at or before the window, end at or after it — the guard production lacks.

Non-spanning reads are not discarded: `getPartialObservationsOfHaplotype`
(`AlleleParser.cpp:3620-3660`) harvests exactly the reads production would silently ref-fill, and
`Samples::assignPartialSupport` (`Sample.cpp:354-420`) decides which candidate alleles each one is
consistent with — testing the bases the read **actually has** as a **prefix** of the candidate's
sequence (the read ran off the right) or a **suffix** of it (ran off the left),
`Sample.cpp:391-395`. One partial read may match several candidates.

**The `1/k` weighting.** `reversePartials` maps one partial read to the *set of candidate alleles it
is consistent with* (`Sample.h:51`, *"for fast scaling of qualities for partial supports"*), so
`k = reversePartials[read].size()` is **how many of the site's candidate alleles this read cannot
tell apart**. The read then contributes `1/k` to each of them — to the depth as well as the
quality, uniformly:

| site | expression |
|---|---|
| `Sample.cpp:41` (`partialObservationCount`) | `scaledPartialCount += 1 / reversePartials[*a].size()` — fractional depth |
| `Sample.cpp:86` (`qsum`) | `qsum += (*a)->quality / reversePartials[*a].size()` |
| `DataLikelihood.cpp:92-93` | `scale = 1/reversePartials[*a].size(); qual *= scale;` |

freebayes' own comment at `DataLikelihood.cpp:97-99`: *"each partial obs is recorded as supporting,
but with observation probability scaled by the number of possible haplotypes it supports."*

**Worked against a six-base window.** Candidates REF `AAAAAA` and a two-base deletion `AAAA`; a read
runs off its own end after three bases, showing `AAA`. That is a prefix of both candidates, so
`k = 2`: freebayes gives **half an observation and half the quality to each**. Production pads `AAA`
with three reference bases, gets `AAAAAA`, matches the REF bucket exactly, and credits **one whole
observation at full quality to REF and nothing to the deletion**.

**Where the two agree is as informative as where they differ.** Take candidates REF `ACGTAC` and a
five-base deletion `A`, with a partial read showing `AC`. `AC` is not a prefix of `A`, so `k = 1` and
freebayes gives it a full REF vote — the same answer production gives, and the right one, because a
read showing `AC` really is evidence against that deletion. So the two callers agree exactly when
the partial read is **informative**, and diverge exactly when it is **not**. That is the worst place
to diverge, and it is the direct consequence of the mechanism: production does not merely skip the
down-weighting, it **manufactures the discriminating bases and then matches on them exactly**, and a
tail copied from the reference can only ever agree with REF.

Complete and partial counts stay separate (`Sample.cpp:11,48`); the channel is on by default
(`Parameters.cpp:435`).

**GATK — no span gate, because it never buckets a read.** Admission is `overlaps`
(`HaplotypeCallerGenotypingEngine.java:348`, `:633`), but the read is scored against each haplotype
by the PairHMM over the bases it actually has. A read ending two bases into a six-base window gets
near-identical likelihoods for REF and DEL, so its likelihood *ratio* is ≈ 1 and it contributes ≈ 0
information — an automatic soft down-weighting, with no counting step to inflate. GATK also carries
the `*` spanning-deletion allele so a site can say "a deletion runs through here" rather than
defaulting to REF.

**bcftools / samtools mpileup — pads the template, never a read.** The indel path realigns each read
over only its real query bases (`bcftools/bam2bcf_indel.c:882-896` clamps `qbeg`/`qend`) and
length-normalises the score (`:547`), so a short read yields a small score *difference* and hence a
small indelQ. Reads entirely inside a deletion get a placeholder large cost (`:912-916`). It does
pad with reference — but on the consensus **template** (`:792-795`), which is safe, never on a
read's observation.

## 5. How it got past the comparison

`doc/devel/reports/reviews/pileup_freebayes_comparison_2026-05-08.md:64-70` lists what the review
deliberately did not examine:

> - Haplotype-window construction (`buildHaplotypeAlleles`, `getCompleteObservationsOfHaplotype`,
>   `getPartialObservationsOfHaplotype`). That is freebayes' *evaluation*-time machinery; we have no
>   analogue and don't need one (Stage 4's grouping plus Stage 5's per-group reconstruction does the
>   same job from `.psp` data).

`getCompleteObservationsOfHaplotype` is not evaluation-time cosmetics — line 3584 *is* the span
gate. The functions were scoped out of the review rather than read.

Separately, the correct behaviour was already written down in this repo and never connected to the
walker — `doc/devel/specs/freebayes_posterior_gt_probs.md:797`: *"when a read only covers part of
the haplotype window, its evidence is distributed across the haplotypes it is consistent with"*.

## 6. No downstream compensation, and the evidence to compensate is discarded

Nothing in `src/var_calling/` treats REF support at a multi-base record differently. `max_ref_span`
([cohort_integration.rs:67-160](../../../../src/var_calling/cohort_integration.rs#L67)) carries span
information but uses it only to align per-sample records onto a common grid, never to modulate
support. Every `num_obs` reader in `per_group_merger.rs` (`:1127`, `:1625`, `:1633`, `:1715`) is a
plain count with no span term, as is the min-observation gate in `variant_caller.rs:378-392`.

And it could not be recovered even if a downstream stage wanted to: `finalise` drops the REF
bucket's chain ids by design ([open_record.rs:158-160](../../../../src/pileup/walker/open_record.rs#L158)),
so the `.psp` carries no per-read identity for REF at all.

## 7. Why it should have the wrong sign for indels — a hypothesis, not a result

The affected reads are precisely the uninformative ones — §4 shows the two callers agreeing whenever
a partial read *can* discriminate and diverging only when it cannot. Production converts a read that
is equally consistent with REF and with a deletion into a **full REF vote**, because the tail it
pads with is the reference by construction. So REF support is inflated at exactly the multi-base
loci where deletions live, which biases against calling them — and `widen` extends the same
inflation to reads that left the walk long before the record reached its final width.

The exposure should scale with **how ambiguous the flanking sequence is**, which is another way of
saying it should be worst in homopolymers and short tandem repeats — where indel placement is
already hardest and where a partial read is least able to discriminate. That is a testable shape,
not just a sign: if the effect is real, the divergent loci should concentrate in low-complexity
sequence rather than spread evenly.

That makes it a candidate mechanism for production's **recorded, still-unexplained indel deficit** —
the one [../../ng/spec/read_preparation.md](../../ng/spec/read_preparation.md) §6 ruled the
soft-mask left-alignment defect out of explaining, since that deficit was measured on GRCh38, which
carries no lowercase bases at all.

**This is a mechanism with the right sign and no measurement behind it.** What would settle it:
stage 2 of the ng port's parity oracle
([../../ng/spec/locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md) §3, §12)
counts, on GIAB HG002 and a tomato CRAM, how many loci, how many reads and how many reference bases
production credits to reads that never sequenced them — and separately the same three for the reads
`widen` extended after they had already expired. If those numbers are small at real depth and read
length, the hypothesis dies cheaply.

## 8. What follows from it

- **Production is frozen; nothing changes there.** This note is the record, and the port-back is
  where a fix would land if the measurement justifies one.
- **ng records partial coverage honestly.** `ReadCoverage` already exists on the shared locus type
  for exactly this ([../../ng/spec/locus_generation.md](../../ng/spec/locus_generation.md) §3), and
  the generic generator computes it from the read's alignment span against the record footprint,
  with a fourth `PartialInterior` variant for the reads a widened record can swallow whole.
- **Whether ng should *consume* partials, and how, is step 7's.** freebayes' prefix/suffix matching
  with `1/|supported alleles|` weighting is the prior art and is MIT-portable. It is the same
  question as the STR path's censored observations
  ([../../ng/spec/locus_generation_ssr.md](../../ng/spec/locus_generation_ssr.md) §8) and should be
  answered once, for both.
- **The `.psp` spec's Step 6 wording should be corrected** whenever production next unfreezes: it
  states an invariant the code does not hold.
