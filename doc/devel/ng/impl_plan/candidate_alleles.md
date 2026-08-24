# ng step 6 — candidate alleles: implementation plan

*Draft, 2026-08-24. Turns the settled design —
[`../spec/candidate_alleles.md`](../spec/candidate_alleles.md) and
[`../arch/candidate_alleles.md`](../arch/candidate_alleles.md) — into build order. **Not a place
for new design:** a question this plan cannot answer from those two goes back to them.*

*This plan builds the shared module and the ordinary SNP/indel path. The repeat-tract path is
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md), which depends on this one and on a field
the merge does not yet carry.*

**Clearing this plan clears the blocker [`calling_loop.md`](calling_loop.md) records twice:** once
selection exists, the calling loop can run on real candidates instead of fixtures, and the GIAB
and tomato regressions the sibling specs name as their definition of done stop being blocked.

---

## Scope

**In:**

- `src/ng/run/cohort_merge/mod.rs` — a `const` path for `MinAltReadShare`, so
  `DEFAULT_ALLELE_SUPPORT` can be a `pub const`; and one doc-comment widening on
  `MinAltReads::reached_by` saying the numerator is the caller's.
- `src/ng/calling/allele_candidates/mod.rs` — `CandidateSelectionConfig` with its two named
  defaults, `SelectionVerdict`, `UnmatchedSupport`, `AlleleRemap`, `LocusSelection`,
  `SelectionScratch`, the private per-allele fold and the ranking comparison.
- `src/ng/calling/allele_candidates/generic.rs` — `select_generic`.
- `examples/ng_candidate_selection_probe.rs` — rewritten to call the module instead of carrying
  its own copy of the rule (Milestone D; it is the measurement oracle).

**Out (later plans):**

- **The repeat-tract path** — [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md). `ssr.rs` is
  not scaffolded here; `mod.rs` declares nothing for it.
- **Nesting `SelectionScratch` inside `CallingScratch`** —
  [`calling_loop.md`](calling_loop.md) A1, which is where `CallingScratch` is built. Here
  `SelectionScratch` stands alone; the loop's A1 makes it a field. Arch §2.4 records why.
- **`GenericSampleEvidence` gaining the dropped-read count, and building the evidence views from
  `remap`** — [`calling_read_likelihoods.md`](calling_read_likelihoods.md) (the type) and
  [`calling_loop.md`](calling_loop.md) E1 (the shaping). Arch §3.2 fixes the hand-off's shape so
  neither has to invent it.
- **Retiring `DiscoveryBar` and routing a discovery round through this rule** —
  [`calling_bakeoffs.md`](calling_bakeoffs.md), which owns `discovery.rs`. Spec §3.4 and §6.3 say
  what it must call and what it must do at the cap.
- **Wiring `select_generic` into the merge's builder for real runs** —
  [`calling_loop.md`](calling_loop.md), which owns that wiring for the loop it sits beside.
- **A bar on summed read quality** — spec Q1's neighbour, spec §3.3. Nothing here reserves a field
  for it; the merge already carries `q_sum`.

## Principles (how the order was chosen)

- **Types first, then implementation** (project rule) — the whole vocabulary in Milestone A, no
  logic until B.
- **The algorithmic heart before the plumbing.** The fold that answers *did this sample support
  this allele* is built and proven on hand-built observations before anything assembles an output
  bundle from it.
- **Reuse over rewrite.** `MinAltReads`, `CandidateAlleles` and the merge's own per-row `q_sum` are
  called as they are; this plan writes a driver, not a second rule (arch §5's table).
- **Isolate the silent steps.** Four steps here produce a *quietly wrong answer* rather than a
  crash, so each lands as **its own commit with its oracle green before and after**, never bundled
  into a neighbour, so a `git bisect` can find it if a number moves: **B1** the per-sample
  denominator, **B2** the ranking, **C1** the remapping, **C3** the leftover sum. Each is marked
  below with the oracle that guards it.
- **Verify against ground truth, and be honest that there is no production oracle.** Production
  has no per-allele bar (spec §9), so byte-parity with it is not available. Two external checks
  replace it: the **measurement oracle** — the standalone implementation in
  `examples/ng_candidate_selection_probe.rs` produced every number the spec quotes, so the shipped
  module must reproduce them (D1) — and the **truth-set property**: over the GIAB trio, of the true
  alternative alleles some sample's reads showed, the bar must keep the fraction spec §3.3
  measured (D2).
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place)

- **The merge, complete** ([`cohort_merge.md`](cohort_merge.md), 16 of 16): `CohortObservation`
  ([`build.rs:922`](../../../../src/ng/run/cohort_merge/build.rs)), `SampleSupport`
  ([`:965`](../../../../src/ng/run/cohort_merge/build.rs)) with its `(allele, read group)` rows
  ([`:1089`](../../../../src/ng/run/cohort_merge/build.rs)) and its partials
  ([`:1130`](../../../../src/ng/run/cohort_merge/build.rs)), `AlleleSupport::q_sum`
  ([`:1245`](../../../../src/ng/run/cohort_merge/build.rs)), and
  `SampleSupport::pooled_support_for` ([`:1070`](../../../../src/ng/run/cohort_merge/build.rs)).
- **The support rule**: `MinAltReads` with `required_of` and `reached_by`
  ([`cohort_merge/mod.rs:424`](../../../../src/ng/run/cohort_merge/mod.rs)).
- **The output table**: `CandidateAlleles` with `new`/`admit`/`bases_of`
  ([`calling/mod.rs:86`](../../../../src/ng/calling/mod.rs)) and `AlleleId`
  ([`types.rs:304`](../../../../src/ng/types.rs)).
- **The measurement harness**: `examples/ng_candidate_selection_probe.rs`, and the benchmark
  inputs it reads — `benchmarks/giab/per_sample/bam/{30x,300x}` with the v4.2.1 truth VCFs, and
  `benchmarks/tomato1/`.
- **Not in place, and not needed here:** `CallingScratch` and `GenericSampleEvidence`. Neither the
  read-likelihoods plan nor the loop plan has started, and this plan depends on neither — selection
  reads the merge's output and writes `CandidateAlleles`, both of which exist.

## Branch and merge

- **Branch** `ng-candidate-alleles`, from `main`. No worktree — the convention for sequential work.
- **Conflict surface, and it is one line.** `src/ng/calling/mod.rs` gains
  `pub mod allele_candidates;`, and [`calling_loop.md`](calling_loop.md) A1 also edits that file.
  Whichever lands second rebases over a one-line addition. Everything else in this plan is new
  files, plus one method and one doc comment in `cohort_merge/mod.rs`, which nothing else in
  flight touches.
- **This plan does not have to wait for the prior or read-likelihood plans** and can run beside
  them.

---

## The steps

### Milestone A — the vocabulary (types, no logic)

**A1. The rule's constants, and a `const` path for the share.**  ✅
`MinAltReadShare` gains a `const` constructor beside its fallible `new`
([`cohort_merge/mod.rs:366-380`](../../../../src/ng/run/cohort_merge/mod.rs)) — its field is
private, so `DEFAULT_ALLELE_SUPPORT` cannot be a `pub const` without one; it panics on a value
outside `0..=1`, which is a compile-time failure for a `const` and so is the right severity.
Then, in the new `calling/allele_candidates/mod.rs`: `CandidateSelectionConfig`,
`DEFAULT_ALLELE_SUPPORT` (floor 2 from `MinAltObs::DEFAULT`, share 5 in 100) and
`DEFAULT_MAX_CANDIDATE_ALLELES = 6`, each with a doc comment carrying its source and marking it
soft. Widen `reached_by`'s doc comment: the numerator is the caller's, and selection and a
discovery round pass different ones. *Depends:* none. *Source:* arch §2.1; spec §3.3, §4.

**A2. The output vocabulary.**  ✅
`SelectionVerdict` (`Selected` / `Truncated { dropped }` / `NotPeriodic`, `#[non_exhaustive]`, the
last documented as repeat-tract-only), `UnmatchedSupport`, `AlleleRemap` with `candidate_for`,
`LocusSelection`, `SelectionScratch`. `calling/mod.rs` gains `pub mod allele_candidates;`. The
two parallelism invariants go in doc comments and are asserted where they are built:
`LocusSelection::unmatched` runs parallel to `CohortObservation::per_sample`, and `AlleleRemap` is
indexed by the merge table's own index. *Depends:* A1. *Source:* arch §2.2, §2.3, §2.4.

> **Checkpoint A:** the module compiles, no logic yet, and every constant names its source and its
> softness. Pause for review.

### Milestone B — the fold (the heart, on hand-built observations)

**B1. The per-sample denominator and the bar. — own commit, do not bundle.**  ☐
The single pass over `CohortObservation::per_sample`: a sample's reads at the locus are the sum of
its rows, **pooled across read groups** and across alleles, which is that sample's compared reads
because the merge admits only complete observations onto alleles; each row then asks
`MinAltReads::reached_by` of its allele. Partials are not read.
**Why this one is isolated:** using depth, or forgetting to pool the read-group rows, gives a
quietly wrong bar at every locus and no test crashes.
*Oracle:* a hand-built locus whose `reads_compared_with_reference` is known independently, carrying
partials, a silent read and one allele shown from two read groups — the denominator must equal the
first and the two group rows must sum rather than the larger winning.
*Depends:* A2. *Source:* arch §3.1; spec §1.3, §3.

**B2. The per-allele summary and the ranking. — own commit, do not bundle.**  ☐
The private `AlleleSummary` (largest within-sample share, samples clearing the bar, cohort read
total, reads and mass for the leftover) filled by B1's pass, and `ranks_above`: share first by
`f64::total_cmp`, then samples clearing, then cohort reads, then the bases.
**Why this one is isolated:** a mis-ordered tie-break or a `partial_cmp` in place of `total_cmp` is
a different truncation at a minority of loci and nothing fails.
*Oracle:* a table built so that each tie-break level in turn is the one that decides, plus the same
table with its rows reversed — the ranking must not move.
*Depends:* B1. *Source:* arch §2.5; spec §4.1.

> **Checkpoint B:** the bar and the ranking are proven on hand-built loci, with nothing assembled
> from them yet. Pause for review.

### Milestone C — the output

**C1. Admission and the remapping. — own commit, do not bundle.**  ☐
`select_generic`'s second pass: seed `CandidateAlleles::new` with the merge's allele 0, `admit`
each surviving alternative in table order, and fill `AlleleRemap` as it goes. A locus where every
alternative failed the bar returns the reference alone with `Selected` — a normal outcome, and the
test says so.
**Why this one is isolated:** an off-by-one in the remapping hands the loop a real but wrong
allele's evidence, which is a wrong genotype rather than a panic.
*Oracle:* the round trip — drop the middle allele of five, then feed every surviving merge row
through `candidate_for` and require the dense id back, and `None` for the dropped one; and the
hand-off of arch §3.2 reproduces the evidence rows exactly.
*Depends:* B2. *Source:* arch §2.3, §3.1; spec §3, §6.2.

**C2. The cap and truncation.**  ☐
Above `max_candidate_alleles`, keep the best by B2's ranking, reference always, and return
`Truncated { dropped }`. Below it, `Selected`. The reference is exempt from both the bar and the
cap. *Depends:* C1. *Source:* arch §2.2, §2.5; spec §4.1.

**C3. The leftover. — own commit, do not bundle.**  ☐
Per sample, sum the merge's own per-row `q_sum` and read count over the alleles that did not
survive — bar or cap — into `LocusSelection::unmatched`. Nothing is re-derived from counts and a
rate. What is *not* in it: partials, reads that produced no observation, reads removed as
evidence, and the reference's reads.
**Why this one is isolated:** a pool that double-counts or misses a row shifts every data
likelihood at that locus, which emission and `QUAL` read as an absolute number, and no genotype
changes.
*Oracle:* the pool equals the bitwise sum of the dropped rows' `q_sum`; a sample with nothing
dropped gets a zero pool with no branch taken to produce it.
*Depends:* C2. *Source:* arch §2.3; spec §5, §5.1.

> **Checkpoint C:** `select_generic` is complete and proven on hand-built loci, including the
> reference-only and truncated outcomes. Pause for review.

### Milestone D — real data, and the two external checks

**D1. The probe calls the module.**  ☐
Rewrite `examples/ng_candidate_selection_probe.rs` to call `select_generic` instead of carrying its
own `summarise`/`keep_by_share` copy, keeping its reporting and its `NG_SELECT_DUMP` /
`NG_SELECT_DUMP_FLOOR` / `NG_SELECT_DUMP_SHARE` surface unchanged. **This is the measurement
oracle:** the standalone implementation produced every number spec §3.3, §4.2 and §5 quote, so the
run must reproduce them — 5,596 alternatives kept of 15,474 on the trio at 300× with the 2% bar,
23 tomato loci of 53,935 above six alleles, 0.36% of tomato's reads in the pool. A difference is a
defect in one of the two implementations and must be traced, not accepted.
*Depends:* C3. *Source:* spec §3.3, §4.2, §5.

**D2. The truth-set property, as a checked-in test.**  ☐
An integration test over a small committed fixture cut from the GIAB trio: of the true alternative
alleles some sample's reads showed, the bar keeps them all but the two spec §3.3 records at 300×.
Written as a property with the count, not as a golden number, so a fixture regenerated at a
different depth still means something.
*Depends:* D1. *Source:* spec §3.3, §12.

> **Checkpoint D:** both external checks green, and the numbers in the spec are reproduced by the
> shipped code rather than by the probe's own copy. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | compiles; every constant carries its source and its softness in a doc comment |
| B | hand-built loci: the denominator against an independently known `reads_compared_with_reference`; two read-group rows summing; each tie-break level deciding in turn; row order reversed with no change |
| C | the remapping round trip with a hole in the middle; the reference-only outcome; `Truncated { dropped }` keeping exactly the top-ranked; the pool as a bitwise sum |
| D | **external**: the probe reproducing the spec's measured numbers on real reads, and the GIAB truth-set property |

**Property tests that outlive the fixtures**, run at every milestone from C on: the surviving list
is a subset of the merge's table, contains the reference, and never exceeds the cap; and adding a
sample that shows only reference reads changes neither the surviving list nor any other sample's
leftover. **The second is spec §3.2's principle as a test — it is the one that fails first if a
cohort term ever creeps into the bar.**

## Out of scope (next plans)

- **The repeat-tract path** — [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md).
- **`SelectionScratch` as a field of `CallingScratch`; the evidence views built from `remap`; the
  wiring into the merge's builder** — [`calling_loop.md`](calling_loop.md) A1, E1 and its own
  out-of-scope list.
- **`GenericSampleEvidence` gaining the dropped-read count** —
  [`calling_read_likelihoods.md`](calling_read_likelihoods.md), with the edit to
  [`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.1 that goes with it.
- **Discovery through this rule, and retiring `DiscoveryBar`** —
  [`calling_bakeoffs.md`](calling_bakeoffs.md).
- **The measurements that would move the two soft numbers** — spec Q1 (the false-allele rate at
  one sample and three reads), Q2 (the cap past 63 samples), Q3 (the share on a second cohort).
  All three are parameters, so no code waits on them.
