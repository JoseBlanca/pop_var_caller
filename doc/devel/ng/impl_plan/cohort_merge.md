# ng cohort merge — implementation plan

**Status:** draft, 2026-08-17. This turns the settled design in
[`../spec/cohort_merge.md`](../spec/cohort_merge.md) and
[`../arch/cohort_merge.md`](../arch/cohort_merge.md) into build order. **It is not a place for new
design** — a step that needs a decision those documents do not make is a gap to take back to them.

The module turns k samples' locus observations into cohort observations, in parallel. It is the
stage that ran on one thread in production and was the wall floor there.

---

## Scope

**In — the direct path only:**

- the three parameters and their defaults;
- the two derivations moved onto the observation types;
- the walk that closes loci by reach and judges them (`close.rs`);
- assembling a survivor: projection onto the locus span, unification into one allele table;
- one builder over one building region (`build.rs`);
- the observation cache and the organiser: overlap resolution, ordered release, eviction
  (`organise.rs`);
- the parallel arrangement, and the proof it changes no output.

**Out:**

- **the psp path's deferred decode** — the module reads whole observations; the two-call source
  arrives with the psp encoding spec ([`../spec/run_streaming.md`](../spec/run_streaming.md) §10).
- **the two caller objects** that drive this module, their look-ahead and their VCF writing —
  [`../arch/run_streaming.md`](../arch/run_streaming.md) §3.
- **calling a cohort observation** — steps 6 to 13, downstream of everything here.
- **where the failed-locus count is reported** — the emission step's
  ([`../spec/cohort_merge.md`](../spec/cohort_merge.md) §13).

---

## Principles (how the order was chosen)

- **The algorithmic heart before the plumbing.** The reach walk and the two verdicts are the
  module; the cache, the organiser and the threads exist to feed and collect them. They are built
  and tested first, against fabricated observations, where no cache can confuse a failure.
- **Serial before parallel, and the serial version is the oracle.** One builder over the whole
  genome produces the cohort observations. Everything the parallel arrangement adds must leave that
  output identical — which is the only assertion that can catch an ownership or overlap defect
  (spec §15).
- **Types first, then implementation**, within each milestone (project rule).
- **Isolate the steps whose failure is silent.** Three produce a plausible wrong answer rather than
  a crash: the verdict order, the allele unification, and the overlap resolution. Each is its own
  commit with its oracle green before and after, so `git bisect` can find it.
- **Reuse over rewrite.** The walk is production's `derive_is_kept` with its columns removed and the
  arithmetic unchanged; the ordered release is `VcfWriter`'s reorder map; the cross-sample merge is
  the argmin shape the read layer already uses.

---

## Preconditions (already in place)

- **The spec and arch are settled**, with one exception named below.
- `SampleLocusObservations` and `SequenceObservation` exist and carry what the module needs —
  `region`, `reference_bases`, `observations`, and the per-read moments
  ([`locus_generation/mod.rs:40,157`](../../../../src/ng/locus_generation/mod.rs)).
- `CensusWriter::add_generic` makes the non-reference comparison inline
  ([`census.rs:2084`](../../../../src/ng/parameter_estimation/joint/census.rs)) — the second caller
  A2 moves onto the shared predicate.
- Production's counterparts are readable and unchanged: `derive_is_kept` and `reach`
  ([`cohort_integration.rs:166-187,46-48`](../../../../src/var_calling/cohort_integration.rs)),
  `PerGroupMerger` ([`per_group_merger.rs:585`](../../../../src/var_calling/per_group_merger.rs)),
  `VcfWriter`'s reorder map
  ([`vcf_writer.rs:162-176`](../../../../src/var_calling/vcf_writer.rs)).

**One decision is still open and it blocks milestone B, not milestone A.** Spec §14 question 2 —
whether a sample's two observations inside one locus combine into one compound allele, and on what
evidence. It must be ruled on before B2. A1 to A4 do not touch it.

---

## The steps

### Milestone A — the walk and the verdicts

Nothing here needs a cache, an organiser, or a thread. Every test fabricates observations directly.

✅ **A1 — the three parameters.** `MaxCohortLocusSpan`, `MinAltObs`,
`CohortLocusBuilderRegionsLen`, each a `NonZeroU32` newtype, with `DEFAULT_MAX_COHORT_LOCUS_SPAN`
= 50, `DEFAULT_MIN_ALT_OBS` = 2, `DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN` = 20.
*Depends:* —. *Source:* [arch](../arch/cohort_merge.md) §1.

☐ **A2 — the two derivations, on the observation types.** `SequenceObservation::matches_reference`,
`SampleLocusObservations::reach` and `::non_reference_reads`; then move
`CensusWriter::add_generic` onto the predicate. **The census's existing tests must stay green** —
that is what says the move changed nothing.
*Depends:* —. *Source:* [arch](../arch/cohort_merge.md) §2.

☐ **A3 — the reach walk, closing loci only.** `LocusCloser` over per-sample cursors: argmin across
their heads, extend the reach, close when the next position falls beyond it. No verdicts yet —
every closed locus comes out. Tests: a SNP at 10 and a SNP at 11 are two loci; a deletion at 10
reaching 14 and a SNP at 12 are one; a late widening pulls in a following observation without a
second pass.
*Depends:* A2. *Source:* [spec](../spec/cohort_merge.md) §4.1; [arch](../arch/cohort_merge.md) §3.

☐ **A4 — the two verdicts, width first.** `Verdict::Failed` / `TooQuiet` / `Build`, decided in that
order. **Own commit, do not bundle.** The order is the silent part: a reference-only chain wider
than `max_cohort_locus_span` must count as failed, not vanish as too quiet, or the count stops
meaning "ground the caller refused". Assert both verdicts on a locus that qualifies for both.
*Depends:* A3. *Source:* [spec](../spec/cohort_merge.md) §3.2, §4.3.

> **Checkpoint A:** loci close correctly and are judged correctly, on fabricated observations, with
> no I/O anywhere. Pause for review.

### Milestone B — assembling a survivor

**B2 needs spec §14 question 2 ruled on.** Do not start it before that.

☐ **B1 — projection onto the locus span.** Each member's sequence widened to the full span, padded
from `reference_bases`. Tests: a narrower SNP inside a deletion's span projects to the span; an
insertion's reference span stays its anchor base.
*Depends:* A4. *Source:* [spec](../spec/cohort_merge.md) §4.2.

☐ **B2 — unification into one allele table.** Identical projections become one allele, the
reference among them. **Own commit, do not bundle.** The silent failure is the opposite of a
crash: one variant written two ways becomes two half-supported alleles, which reads as a noisy
site. Test a deletion presented at two placements — it must unify, and the test must state that it
relies on left-alignment upstream.
*Depends:* B1, and the §14 Q2 ruling. *Source:* [spec](../spec/cohort_merge.md) §4.2.

☐ **B3 — `CohortObservation` and `SampleSupport`.** Per sample, support against the allele table,
moments summed where two of its own observations projected onto the same allele; support never
merged across alleles. A sample with no coverage has no support, which stays distinct from
reference-only support.
*Depends:* B2. *Source:* [arch](../arch/cohort_merge.md) §4.

> **Checkpoint B:** a cohort observation is built from fabricated per-sample observations, with one
> allele table and per-sample support. Pause for review.

### Milestone C — one builder, serially, over the whole genome

☐ **C1 — `build_region`.** Walk, judge, assemble; return the survivors and the failed spans.
*Depends:* B3. *Source:* [arch](../arch/cohort_merge.md) §4.

☐ **C2 — a serial driver, and it becomes the oracle.** One builder, no cache, no organiser, reading
observations straight from a source over the whole analysed stretch. Its output is what every later
milestone must reproduce.
*Depends:* C1. *Source:* [spec](../spec/cohort_merge.md) §15's partition-invariance oracle.

> **Checkpoint C:** cohort observations from real per-sample observations on a small fixture, one
> builder, single-threaded. Pause for review.

### Milestone D — the observation cache

☐ **D1 — `ObservationCache`.** One forward reader per sample; `cover`, `with_observations`,
`evict_before`. Tests: a window covering a region also holds an observation that started before it
and reaches in; eviction drops nothing a live region can still reach.
*Depends:* C2. *Source:* [spec](../spec/cohort_merge.md) §6.4; [arch](../arch/cohort_merge.md) §4.

☐ **D2 — the serial driver reads through the cache.** Output byte-identical to C2.
*Depends:* D1. *Source:* [spec](../spec/cohort_merge.md) §6.4.

> **Checkpoint D:** the cache feeds the serial builder and changes nothing. Pause for review.

### Milestone E — the organiser and the parallel arrangement

☐ **E1 — ordered release, no overlaps yet.** `Organiser::submit` / `drain_ready` / `is_finished`,
keyed by region index, drained on next-expected; exactly one outcome per region including empty
ones; a gap is an error rather than a truncation.
*Depends:* D2. *Source:* [spec](../spec/cohort_merge.md) §6.3; [arch](../arch/cohort_merge.md) §4.

☐ **E2 — overlap resolution.** Of two overlapping loci the earlier start stands; a failed locus's
span displaces the same way an emitted locus does. **Own commit, do not bundle.** The silent
failure is a locus built from a partial picture surviving into the output because nothing displaced
it — which needs a wide deletion beginning before a building region and reaching into it, so the
fixture must contain one.
*Depends:* E1. *Source:* [spec](../spec/cohort_merge.md) §6.1.

☐ **E3 — builders in parallel.** Regions handed out, builders reading the shared cache, the
organiser collecting.
*Depends:* E2. *Source:* [spec](../spec/cohort_merge.md) §6.1, §6.2.

☐ **E4 — assert the milestone.** Same cohort observations and same failed spans as C2, at 1, 2, 4,
8 and 16 builders and at several values of `cohort_locus_builder_regions_len`.
*Depends:* E3. *Source:* [spec](../spec/cohort_merge.md) §15.

> **Checkpoint E:** parallel output identical to the serial oracle, at every builder count and
> every region width. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | fabricated observations: the grouping cases of spec §4.1, and a locus that qualifies for both verdicts taking the width one |
| B | a deletion written at two placements unifying to one allele; per-sample support summed within an allele and never across |
| C | cohort observations on a real fixture, single-threaded — the oracle every later milestone reproduces |
| D | output byte-identical to C2 with the cache interposed |
| E | output identical to C2 at 1/2/4/8/16 builders and several region widths, including a fixture with a deletion reaching across a building-region boundary |

**The north-star oracle is C2**, and it is deliberately the simplest thing that produces the right
answer: everything after it is about speed and memory, so anything that changes the answer is a
defect rather than a trade.

---

## Out of scope (next plans)

- **The psp path's two-call source and its deferred decode** — the psp encoding spec
  ([`../spec/run_streaming.md`](../spec/run_streaming.md) §10). Nothing in this module changes when
  it lands: the walk consumes a position, a reach and a count per observation either way.
- **The caller objects, the look-ahead and the VCF writer** —
  [`../arch/run_streaming.md`](../arch/run_streaming.md) §3, which owns what drives this module.
- **Calling a cohort observation** — steps 6 to 13.
- **The measurements that set the two starting values** — `cohort_locus_builder_regions_len` at 20
  and `max_cohort_locus_span` at 50, spec §14 questions 1 and 3. Both are parameters, so the code
  does not wait on them; the sweep wants E3 to exist first.
