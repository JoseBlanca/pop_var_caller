# ng parameter pre-pass, the STR path (step 4) — implementation plan

**Status:** draft, 2026-08-07, revised 2026-08-09 for spec §4.5 (a second floor on slipped reads,
the split provenance it forces, and two new anchors at Milestone G); **revised 2026-08-11** for the
loci's new provenance and one new milestone. Three things moved and none of them is the
mathematics: the loci reach this step from a repeat catalog built once per reference rather than
from a scan run while the sample is read (spec §5.3), so every step that touches real data takes a
catalog as an input; the copy floors became `[8, 6, 6, 6, 5, 4]` on 2026-08-10, which re-routed nine
in ten non-mononucleotide loci to the generic path and dated every measured number this plan used to
assert (spec §5.1.1); and **Milestone H asks what a locus that is not the tract its stratum thinks
it is does to the four numbers** (spec §4.6), which is the question the sibling path had to answer
after its own end-to-end run and answered with a second class of site. The build order for the **STR
half of step 4**: one accumulator keyed
by `(read group, motif period, reference repeat count)`, the four numbers fitted from it — how often
a read slips, which way, how far, and a per-base substitution rate — and the summary a person reads
instead of several hundred per-stratum records. Design is settled in
[`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) (spec) and
[`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) (types & interfaces), on the
shared framing of [`parameter_prepass.md`](../spec/parameter_prepass.md) and under the shared arch
docs ([step interfaces](../arch/ng_step_interfaces.md), [module layout](../arch/module_layout.md)).
This turns that design into build order; it is **not** a place for new design.

**It follows [`parameter_prepass_generic.md`](parameter_prepass_generic.md)**, which builds the
`fitting/` module this plan is the second consumer of. That plan named the two changes this one asks
of `fitting/`, and neither is a rewrite: `fit_mixture_weights` widens past three genotypes, and a
multi-start maximiser lands beside `fit_by_profile_scan` rather than replacing it. **Whether the seam
was cut in the right place is what this plan finds out**, and Milestone D is where.

## What this step can be checked against

The sibling plan's answer holds here and needs no restating: there is no production oracle, because
production's STR pre-pass pools reads from loci that passed a confident-genotype gate
([`src/ssr/cohort/prepass.rs`](../../../../src/ssr/cohort/prepass.rs)) — the exact bias this step
exists to remove, so agreeing with it would be the bug. There is no legible fixture either: every
number here is the argmax of a sum over hundreds of thousands of loci.

**So the oracle is the estimator's bias, computed exactly** — replace each cell's observed count with
its probability under a known truth, maximise, and the answer is what an infinite genome returns,
with no sampling noise in it. **The harness that does this already exists and is green**:
[`examples/ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs), written up in
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§6. Every milestone below is proven against it, against an algebraic identity, or against HG002's
truth genotypes — never against itself.

**Three cheap checks come before any fit, and each rejects a broken scoring rule in one line** (spec
§10.1): the rule sums to one over the entry space at any parameter values; no bucket is charged a
negative number of reads; and with the slippage level at zero every locus's reads land on its own
alleles. **Any change to the scoring rule or the entry key re-runs all three first.**

**One diagnostic that looks sufficient and is not, because it has already produced a retracted
finding.** A search that returns the same answer from four different starting points reads as a
search that found something. It is not: a deterministic optimiser returns the same point from every
start wherever the objective is flat, which is exactly what happened on 2026-08-06 — four starts
disagreeing tenfold about the fall-off all returned it to four decimal places, and it was the inner
climb stopping short (research note §6.3.1). **Pair every spread check with the score at the truth**,
which costs one evaluation and which a correctly specified model cannot be beaten at. A fitted point
scoring *above* the truth is a defect in the test, not a finding about the estimator.

---

## Scope

**In:** `src/ng/parameter_estimation/ssr/` — `mod.rs`, `locus_offsets.rs`, `stratum_table.rs`,
`slippage.rs`; the STR vocabulary added to `types.rs` (`SsrPeriod`) and the step-4-local scalars; the
sparse table of locus shapes with its merge; the slippage noise model and its marginal end-bucket
rule; the substitution rate's closed form; the multi-start search; borrowing against **both** floors
— loci for the whole model, slipped reads for the direction split and the fall-off alone (spec
§4.5) — the monotonicity walk, and the per-read-group summary; both entry points. Plus the two
additive changes to `fitting/`, and the diagnostic that reports how much of a stratum sits in shapes
its own fit cannot explain, with the measurement that says what such loci cost (Milestone H).
**One test lands outside this folder**: the agreement between this path's substitution rate and the
generic path's belongs at step 4's own surface, which is the only place holding both (G5).

**Out (each handed to a named owner):**

- **The generic half of step 4** — [`parameter_prepass_generic.md`](parameter_prepass_generic.md),
  which this plan depends on rather than duplicates. Its Milestones A and D are preconditions below.
- **Where the mismatch count for the substitution rate comes from** — `OPEN` in arch §2.3, and the
  better of the two answers touches step 3's aligner. **This plan builds the answer that needs
  nothing upstream** (compare each read against the motif tiled to the length that read shows) and
  records its caveat in the emitted provenance. Adopting the aligner's own count later changes one
  function's body and no signature. **Home:** whichever plan revises step 3's alignment output.
- **The per-period copy floors** — settled, and settled *against* this plan's measurements rather
  than by them: `[8, 6, 6, 6, 5, 4]` since 2026-08-10 (spec §5.1.1,
  [`segment_criteria.rs:444`](../../../../src/ng/region_typing/segment_criteria.rs)). They decide
  which tracts arrive here, not what is done with them. **Nothing in this plan waits on them, and
  every measured number it quotes predates them** — which is why the real-data steps assert shapes
  and record numbers rather than the other way round.
- **The STR census, and the per-locus model it would make askable** — spec §8.5, §8.6; the census
  has a spec ([`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md)) and
  no architecture document.
- **Reconciling this path's per-stratum allele spectra with the cohort's STR diversity** — spec §6.
  **Home:** [`parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md) §3.
- **Long-allele recovery and the library profiling it needs** — spec §7, deferred on coverage.
- **How many samples are walked at once** — a driver decision, and the same one the sibling plan
  defers.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The mathematics before the plumbing.** The scoring rule and the search (Milestone D) are built
  and proven against the harness's exact answers **before a single locus is read**. What must be
  right is the model; the accumulator only has to hand it the right counts.
- **A directly-filled table is the oracle for the accumulator that fills it.** Milestone D's fits run
  on tables built entry by entry from known parameters. Only once they recover the truth does a
  locus stream fill one (F3). An accumulator bug and a model bug then cannot hide each other.
- **Isolate the steps whose failure is silent, and say so.** Most of this module fails loudly — a
  panic, a failing test. **Six do not**: the complete-witness rule (C2), the read cap's draw (C3),
  the end-bucket scoring (D2), the multi-start search (D4), the monotonicity merge (E4) and the
  contamination measurement (H2). Each returns a plausible number nobody can check, and **four of
  them have already produced a wrong one** during the measurement work. They land as **their own commit with their oracle green before
  and after**, so `git bisect` can find one if a parameter later moves. They are marked **own commit,
  do not bundle**; no other step is.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); a native host build at
  completion.

## Preconditions (already in place, or named as the gate)

- **The generic plan's Milestones A and D are done.** This unit consumes `ErrorRate`, `Ploidy`,
  `PloidyMap`, `Estimate<T>` and `Provenance` (generic arch §2.1, §2.4, §3), the `NoiseModel` trait
  and `fit_mixture_weights` (generic arch §4). **It does not need Milestones B, C, E, F or G** — the
  generic cell table, its accumulators, its four fits and its anchors are all beside the point here.
  So this plan can start once the sibling reaches its Checkpoint D.
- **The STR locus stream exists and runs on real alignments.** `SsrGenerator`
  ([`locus_generation/ssr.rs`](../../../../src/ng/locus_generation/ssr.rs)) yields
  `LocusKind::Ssr` loci with `complete_observations()`, `reference_bases` and per-observation
  `read_group` — the only things this step reads from a locus (arch §6). Two example programs already
  drive it over the tomato and HG002 alignments end to end
  ([`ng_str_table_memory.rs`](../../../../examples/ng_str_table_memory.rs),
  [`ng_str_stutter_by_library.rs`](../../../../examples/ng_str_stutter_by_library.rs)).
- **A repeat catalog exists for each reference the real-data steps use, and it is an input like the
  reference is.** Nothing finds a tandem repeat while a sample is read any more (spec §5.3): the
  regions the locus stream walks are read from a file written once by
  `pop_var_caller_exp repeat-catalog`, and a driver of this step opens it. **So Milestones F and G
  do not run without a catalog for tomato SL4.0 and for GRCh38**, built at the catalog's own floors
  (`[5, 5, 4, 4, 4, 3]` copies, 15 bp flank, periods 1 to 6), which sit under every floor a run
  routes on. A catalog is bound to the reference it was built from and refuses to be paired with
  another. Both files exist on the development machine under `tmp/`; **that is scratch, so a fresh
  machine builds them first**, and neither the build nor its cost belongs to this plan.
- **The numbers this plan used to assert are dated, and the milestones below say so where they use
  them.** The copy floors moved on 2026-08-10 (spec §5.1.1), so a walk today admits a smaller
  population of non-mononucleotide loci than every measurement quoted here was taken over. Where a
  step asserted a measured count, it now asserts the shape and **records the current number when the
  real accumulator first walks each cohort**.
- **The exact-bias harness runs and its numbers are recorded.**
  [`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs), research note §6. It
  is the oracle for Milestones D and E and must be green before either starts.
- **Most of Milestones B and D already exist as working code, and this plan ports rather than
  re-derives them.** The library survey built them to answer the copy-floor question:
  [`examples/shared/stutter_model.rs`](../../../../examples/shared/stutter_model.rs) carries the slip
  kernel with its truncation renormalised, the per-genotype bucket probabilities with the marginal
  end-bucket rule **and** the edge plug-in beside it as a control, the climb over the genotype
  frequencies at arbitrary genotype counts, the four starting points and the multi-start driver, and
  the three algebraic gates; [`examples/shared/stutter_table.rs`](../../../../examples/shared/stutter_table.rs)
  carries `LocusShape`, the sparse table, `not_whole_repeat_share`, the allele support with its
  low-end clip, and the per-stratum fit. **Port them; do not write the same mathematics twice.** The
  reason is not economy: four of this plan's five silent-failure steps produced a wrong number once
  already, and the end-bucket rule is where two of them were, so a second independent derivation is a
  second chance at the same mistake rather than a check on it. **What a port does not inherit is the
  acceptance gate** — every oracle below still has to pass against the ported code, and D2's four
  identity checks are what prove the port faithful.
- **No production dependency.** `src/ssr/` and `src/pileup/` are frozen. Production's stutter
  pre-pass is **not** reused (arch §6), for the reason under "What this step can be checked against".

---

## The steps

### Milestone A — vocabulary and the local types (types, no logic)

**A1. Scaffold `parameter_estimation/ssr/`.**  ✅
`mod.rs`, `locus_offsets.rs`, `stratum_table.rs`, `slippage.rs`, each with its `#[cfg(test)]` block;
wire `pub mod ssr;` into `parameter_estimation/mod.rs`. A folder rather than a file for the reason
`generic/` is one: the shaping of data and the mathematics on it never share a file. **No trait over
the accumulator** — nothing generic drives it, and the walk knows which object it is filling.
*Depends:* generic A1. *Source:* arch §Module home, [module layout](../arch/module_layout.md).

**A2. `SsrPeriod` into `types.rs`; `RepeatCount` and `Stratum` beside it.**  ✅
`SsrPeriod` is a checked `u8` rejecting zero and anything above `MAX_MOTIF_LEN`, because a period of
zero divides by zero when a tract length becomes a repeat count. It is shared vocabulary with steps 6
and 7, which is why it goes in `types.rs`; `Motif` gains `ssr_period()` **beside** its existing
`period()` rather than changing that accessor's type under its callers
([`types.rs:369`](../../../../src/ng/types.rs)). `Stratum { period, repeats }` orders by
`(period, repeats)` so the monotonicity walk visits neighbours in order. Unit tests: period 0 and
period 7 rejected; `Stratum` ordering is by period then repeat count. *Depends:* A1.
*Source:* arch §2.1.

**A3. The offset scalars and the two widths.**  ✅
`WholeRepeatOffset`, `OffsetBucket`, `OFFSET_HALF_RANGE = 4`, `OFFSET_BUCKETS = 9`,
`ALLELE_OFFSET_LIMIT = 6`, `MAX_LOCUS_READS = 12`, `GUARD_SHARE_LIMIT = 0.10`. **The two widths are
different things and only one is load-bearing**, which the constants' doc comments must say: the
*recorded* offset range can be narrow, because an end bucket scored by its marginal still returns the
slippage level to within 0.05% at ±1 against alleles reaching ±3; the *allele support* the fit may
place mass on is what decides the answer, and it clips at the low end so a stratum at 4 repeats
reaches only −4. Unit tests: `bucket_of` is total and monotone over the offset range and saturates at
both ends; the allele support at repeat count 3, 6 and 20 has 10, 13 and 13 lengths.
*Depends:* A2. *Source:* arch §2.1.

**A4. The three slippage rates and the model that holds them.**  ✅
`SlipRate`, `SlipGainShare`, `SlipStepDecay` — three types and not one shared `Probability`, because
they are all fractions in `[0, 1]` and one type would let a direction split be handed to something
expecting a slippage rate and compile. Each copies `MismatchFraction`'s shape
([`types.rs:243`](../../../../src/ng/types.rs)): private field, `try_new`, `.get()`. `SlippageModel`
holds all three. Extend `DomainError` with their three variants — its doc already says later
constrained types add their own ([`types.rs:268`](../../../../src/ng/types.rs)). Unit tests:
boundaries accepted, out-of-range rejected. *Depends:* A2. *Source:* arch §2.1, §2.4.

**A5. The output types.**  ✅
`StratumFit`, `SlippageStart`, `SsrSampleParameters`, `StratumFitSummary`, `SsrAccumulationCounts`.
Types only. `StratumFit` carries **two** provenance lists and not one — `fitted_over` for the level,
`shares_fitted_over` for the direction split and the fall-off — because those starve at rates 20,000
apart and a stratum routinely measures the first well while borrowing the second (spec §4.5); plus
`slipped_reads`, without which a level of 0.0003 standing on 4 reads is indistinguishable downstream
from one standing on 4,000. A borrowed or merged value is a different claim from one fitted in place
and a consumer must be able to tell. *Depends:* A4. *Source:* arch §2.4, §4.3.

**A6. `SsrEstimationError`.**  ✅
`NoFittableStratumAtPeriod`, `SlippageNotIdentified`, `Domain`, with `MIN_LOCI_TO_FIT = 1_000`,
`MIN_SLIPPED_READS_TO_FIT_SHARES = 4_000` and `START_AGREEMENT_LIMIT = 1.06`. **The second floor's
doc comment must carry its derivation**, because it looks arbitrary and is not: at §3's measured
values, holding the direction split to 6% of itself takes about 1,400 slipped reads and holding the
fall-off to the same takes about 4,000, so the fall-off binds. **This path's own enum, not variants bolted onto the generic path's**:
the two units fail differently — a thin stratum has neighbours to borrow from and a sample's
heterozygosity does not. `NoFittableStratumAtPeriod` deliberately has **no default value to fall back
on**, because a slippage rate spans twenty-two-fold across repeat counts within one dataset, so any
constant would be wrong for most strata. Unit test: each message names the sample and the number that
was too small. *Depends:* A4. *Source:* arch §4.2.

> **Checkpoint A:** the vocabulary compiles; every constrained rate rejects what it must; the two
> widths are pinned by test, including the low-end clip. Pause for review.

### Milestone B — the table of locus shapes (storage, no loci)

**B1. `LocusShape` and its invariant.**  ✅
Nine bucket counts and a guard count, all `u8`, with `counts.iter().sum() + not_whole_repeat ==
depth` holding always. Ordered and hashable, so it can key a table and so iteration order is fixed —
which is the whole of the determinism requirement. **The depth is exact, not binned**, and that
follows from the cap rather than from a separate decision: `MAX_LOCUS_READS` bounds it at a dozen
values, so there is nothing for a ladder to save. Unit test: a shape whose counts exceed the cap
cannot be built. *Depends:* A3. *Source:* arch §2.2.

**B2. `StratumTable` — storage, `add_locus`, `shapes`, `loci`.**  ✅
A `BTreeMap<LocusShape, u32>` and two `u64` composition counters. **Sparse and not dense, and that is
forced rather than chosen**: an entry is a whole locus's split across ten buckets, so the possible
space is 220 shapes at three reads a locus and 293,930 at twelve, of which only a small
data-dependent corner is ever occupied. `BTreeMap` and not `HashMap`, because every fit is a
floating-point sum over entries and floating-point addition is not associative. `shapes()` returns a
`Vec` and not an iterator, because the search re-walks the entries once per candidate. Unit tests: two
loci with the same shape make one entry with a count of two; `shapes()` is stable in order across
runs. *Depends:* B1. *Source:* arch §2.2.

**B3. `merge`, `substitution_rate` and `not_whole_repeat_share`.**  ✅
`merge` is element-wise integer addition, so it is associative and exact and shards merge to the table
of the union. `substitution_rate` is mismatched over compared — **a division, not a search** (spec
§4.1), and the closed form is the maximum rather than a moment estimate. `not_whole_repeat_share`'s
denominator is the reads differing from **the reference tract length**, not from the allele, which
the accumulator cannot know; the share is therefore diluted relative to the model's and never
inflated, so a stratum crossing the limit on it has crossed it on the model's too. Unit tests: a
table split arbitrarily in two and merged equals the unsplit one, entry for entry, in either merge
order; the substitution rate of a table with no bases compared is `None` rather than zero.
*Depends:* B2. *Source:* arch §2.2, spec §4.1.

> **Checkpoint B:** the table stores, merges exactly and reports its two diagnostics; the
> substitution rate is proven to be the maximum and not merely a ratio. Pause for review.

### Milestone C — one locus → one entry (data shaping)

**C1. `stratum_of`.**  ✅
The reference tract's period and repeat count, from `reference_bases.len() / motif.period()` — both
of which the locus carries. **A pure function of the reference**, which is what makes every sample
stratify identically so a cohort can compare strata. A tract whose reference length is not a whole
number of copies is **counted and skipped, not rounded**. Unit tests over hand-built loci: a
non-`Ssr` locus returns `None`; a 13-base tract at period 3 is counted and skipped.
*Depends:* B1. *Source:* arch §2.3.

**C2. `shape_of` — complete witnesses only.**  ✅ **Own commit, do not bundle.**
A read's offset is `(observation.bases.len() − reference_bases.len()) / period`, whole only when that
difference divides by the period; otherwise the read goes to the guard bucket. **The silent failure
this isolates:** a partial witness saw only part of the tract, so its length is a **lower bound** —
scoring it as a length reads as a read that lost repeats, which is a direct bias in the direction
split, the one parameter §3 of the spec exists to protect and the one that inverted on real data.
`reads_without_observation` does not enter the depth; those reads covered the tract and witnessed
nothing. `reads_discarded_by_cap` does **not** skip the locus — the generator's own reservoir is a
random subsample, so a locus it fired on is a locus observed at a lower depth — but it is counted,
because a run where it fires everywhere is a run whose depths are the cap's and not the data's.
*Oracle:* a hand-built locus whose partial witnesses all show short tracts must produce a shape with
every read at the origin, and the same locus scored without the guard must produce a visibly
different one — so the test is proven to bite. *Depends:* C1. *Source:* arch §2.3, spec §3.

**C3. The read cap — subsample, seeded from the locus's position.**  ✅ **Own commit, do not bundle.**
A locus deeper than `MAX_LOCUS_READS` is entered from a uniform random subsample of its reads down to
the cap. A subsample is exact rather than approximate: thinning a locus's reads uniformly leaves the
bucket counts distributed exactly as they would be at the lower depth. **The silent failure this
isolates:** the seed. Seeded from the locus's position, a region-sharded walk and a single-threaded
one keep the same reads and `merge` stays an equality; seeded from anything else — a counter, the
thread, the clock — they diverge, and the divergence is a few reads per deep locus, which no test
that does not compare two walks would ever show. *Oracle:* the same locus gives the same draw on
every run and in every shard layout, and over many positions the kept bucket counts are
hypergeometric in mean and variance. *Depends:* C2. *Source:* arch §2.1, spec §4.1.

**C4. `composition_of` — bases compared and bases mismatched.**  ✅
Compare each read's bases against the motif tiled to **the length that read shows**, so a mismatch is
a substitution and not a slip. **This is an alignment, not a call**, so it does not reintroduce the
threshold-then-count bias step 4 exists to remove. It is the answer that needs nothing upstream; the
better one — having the aligner emit the count it already computes while scoring
([`alignment/emission.rs:43`](../../../../src/ng/alignment/emission.rs)) — touches step 3 and is out
of scope above. **Its caveat travels in the emitted provenance**: an impure tract's interruptions are
charged to the substitution rate, which `SsrSegment::purity_fraction()`
([`segment_criteria.rs:209`](../../../../src/ng/region_typing/segment_criteria.rs)) makes measurable
per stratum. Unit tests: a perfect tract mismatches nothing; a tract with one interior substitution
mismatches once at every length the read shows. *Depends:* C1. *Source:* arch §2.3.

**C5. `SsrAccumulators`, `add_locus`, `merge`, `adjustments`.**  ✅
One `StratumTable` per `(read group, stratum)`. `add_locus` **borrows** the locus and passes it on
untouched, ignores a `kind` that is not `LocusKind::Ssr`, and tallies rather than repairs. **A locus
covered by two read groups makes two entries and that is sound** — the genotype is drawn once for the
locus and enters both through the same mixture, so the product over them is a composite likelihood
and the split costs precision, not correctness. What must not be split is a locus's reads *within*
one read group, which the entry key prevents. `loci_without_whole_repeat_reference` **must read near
zero** and is a bug report against the classification the catalog reader applies if it does not — a
tract admitted as a locus whose reference length is not a whole number of copies is a delimiting
fault upstream, never something this unit rounds away. Unit tests: a non-STR locus changes
nothing; three shards merged in every order give identical tables and identical counters.
*Depends:* C3, C4, B3. *Source:* arch §4.

> **Checkpoint C:** loci reduce to entries, the cap is reproducible from position alone, and sharded
> accumulation is proven order-independent as an equality rather than a tolerance. Pause for review.

### Milestone D — the noise model and the search (the mathematics, no loci)

**D1. The slip kernel.**  ✅
`P(a read shows exactly d whole repeats more than its allele)` from the three slippage parameters: no
slip with probability `1 − level`; otherwise a direction drawn from the gain share and a distance
drawn from a geometric fall-off, **one fall-off shared by both directions** (spec §3). The truncation
at the largest representable step is renormalised so it loses no mass. Unit tests: the kernel sums to
one over its support at every parameter setting tried; at a level of zero it is one at `d = 0` and
zero elsewhere; the ratio of the two-step to the one-step term equals the fall-off in both
directions. *Depends:* A4. *Source:* spec §3, harness `Slip::p`.

**D2. `SsrNoiseModel::genotype_likelihoods` — the marginal end-bucket rule.**  ✅ **Own commit, do
not bundle.**
A genotype is an unordered pair of allele offsets. Each read picks one of the two copies and then
slips, so a bucket's probability is the average of the two copies' slip kernels, and a shape's
probability is one multinomial over the buckets. **An end bucket's probability is the sum over every
offset it absorbs, never the probability of sitting exactly on the edge.** **The silent failure this
isolates:** the plug-in is the tempting shortcut and it is wrong twice over — scoring at the edge
costs **52% of the slippage level** where the alleles reach ±3 and the recorded range is ±1, and
rescaling it to sum to one protects the fall-off while costing **+33% of the level** where 30 in 100
slipped reads take a second step, which is the regime long tracts sit in. Neither shows on the
outside. *Oracle:* the three algebraic gates first, each one line and none needing a fit — the rule
sums to one over the entry space at any parameters (the un-rescaled plug-in sums to 0.9488 at ±1 and
is rejected here, before any fitting); no bucket is charged a negative number of reads; and at a
level of zero every locus's reads land on its own alleles. Then agreement with the harness's
`genotype_bucket_probs` to floating point on one world's entry space. **Any later change to this
expression re-runs all four.** *Depends:* D1, generic D2's `NoiseModel` trait. *Source:* arch §3,
spec §4.1, research note §6.4.

**D3. Widen `fit_mixture_weights` past three genotypes.**  ✅
Its declared return type is `SmallVec<[f64; 3]>`, which is the diploid generic path's genotype count.
A stratum here has up to 91 — thirteen allele lengths at `ALLELE_OFFSET_LIMIT = 6`, fewer where the
support clips at the low end — so the return type widens to a `Vec<f64>` or becomes generic in its
inline capacity. **One line, in the shared module rather than copied here**, and the first of the two
changes the sibling plan anticipated. The climb itself is unchanged: it is the same concave problem
and the same code, and convergence failure stays a bug rather than a data condition. Unit test: the
existing generic three-genotype test still passes, and a hand-built 45-genotype table recovers its
known weights from any interior start. *Depends:* generic D1. *Source:* arch §3.

**D4. `fit_by_multistart`, and the spread it must report.**  ☐ **Own commit, do not bundle.**
Maximise the three slippage parameters from several starting points, climbing the genotype
frequencies at every trial, and return the best-scoring **with every start's outcome beside it**.
`SLIPPAGE_STARTS` is four starts that disagree about the level, the direction and the decay **at
once** — the level as a multiplier on a moment estimate rather than an absolute value, because a
stratum's rate spans twenty-two-fold across repeat counts and a fixed ladder of absolute rates would
start every stratum in the wrong place. **The silent failure this isolates, and it is the one that
cost a published finding:** starts that agree are not evidence of a fit. A deterministic search
returns the same point from every start wherever the objective is flat, so four starts agreeing to
four decimal places is exactly what a search that never looked also produces. It is also the trap the
generic path's inbreeding fit fell into from the other side — five starts disagreeing about the
headline number while sharing one guess at a nuisance axis returned a confident zero on a genome 29%
covered by runs. *Oracle:* the harness's control, and it must read **exactly** zero — generate and
fit under the same key with the reference origin, and get 0.000% on the level, 0.0000 on both shares,
four starts agreeing to 1.000×. Paired with the score at the truth, which the fitted point may not
exceed. **Convergence failure here is a data condition, unlike the inner climb**: the outer search
has no concavity proof, so it is capped, the best-scoring iterate kept, and the termination reported.
*Depends:* D2, D3. *Source:* arch §3, spec §4.2, research note §6.2, §6.3.1.

> **Checkpoint D:** the scoring rule passes all four identity checks, the search recovers a known
> truth to exactly zero bias, and the score at the truth is unbeaten. **Nothing has read a locus
> yet.** This is also where the `fitting/` seam is judged: if D2 and D4 needed more from it than the
> two changes above, say so here rather than working around it. Pause for review.

### Milestone E — the four fits, in order

**E1. The substitution rate.**  ☐
Mismatched bases over bases compared, per stratum. One division, and it needs none of the other
three, which is why it goes first. Where a stratum holds reads of two different true rates the pooled
counters return their base-weighted mean, which is the right answer for a model carrying one rate.
*Oracle:* the harness recovers 0.0030 from a truth of 0.0030 by a search that had no need to run.
*Depends:* B3, C5. *Source:* arch §4.1, spec §4.1.

**E2. The three slippage parameters, per stratum.**  ☐
`fit_by_multistart` over that stratum's shapes, with `fit_mixture_weights` climbing the genotype
frequencies at each trial. Genotype frequencies are **fitted freely** over unordered allele pairs
rather than tied through one allele frequency, matching the generic path and for the same reason: a
Hardy-Weinberg tie presumes the inbreeding coefficient is zero, and the inbreeding coefficient is a
quantity this run measures rather than assumes. Each fit records its starts and their scores in
`StratumFit::starts_tried`, and raises `SlippageNotIdentified` when they span more than
`START_AGREEMENT_LIMIT` in the level. *Depends:* D4, E1. *Source:* arch §4.1, spec §4.2.

**E3. Borrowing for a thin stratum — two floors, and the fit splits between them.**  ☐
Below `MIN_LOCI_TO_FIT`, take the neighbouring repeat counts at the same period rather than fitting
noise, marked `Provenance::Borrowed` with `fitted_over` naming the strata it came from. A period with
no fittable stratum anywhere raises `NoFittableStratumAtPeriod` rather than defaulting.

**Then the second floor, which is on slipped reads and not on loci** (spec §4.5). A stratum that
clears `MIN_LOCI_TO_FIT` can still put single digits behind the direction split and the fall-off,
because those are measured only by the reads that moved: at the 0.091% level §5 measures below four
repeats, 100,000 loci at 5 reads each yield 455 slipped reads, 77 gains and about 5 second-step
gains. **The level is untouched by the same thinness** — a proportion over all 500,000 reads, good
to about 5% of itself — so below `MIN_SLIPPED_READS_TO_FIT_SHARES` the level is kept and only the
two shares are borrowed, with `shares_fitted_over` naming their source. Emit `slipped_reads` beside
the fit either way. Unit tests: a thin stratum between two thick ones borrows the whole model and
says so; a stratum with a million loci and 40 slipped reads keeps its level and borrows its shares,
with the two provenance lists differing; a period whose every stratum is thin errors.
*Depends:* E2. *Source:* arch §4.1, §4.2, spec §4.5.

**E4. The monotonicity walk — merge and refit.**  ☐ **Own commit, do not bundle.**
Last, because it reads the fitted sequence. Visit each period's strata in repeat-count order; where a
fitted level falls below its predecessor's, merge the two tables and refit, repeating until the
sequence rises. **The silent failure this isolates:** a merge **changes the estimate** and does so
without failing anything. Two strata pooled return close to the loci-weighted mean of their levels,
so each then carries its own distance from it — a 1.5-fold difference between neighbours costs about
a quarter of the level, a two-fold difference about half, a four-fold difference up to 141%. On real
strata slippage rises about 1.3-fold per repeat count, so one merge costs on the order of 15 to 25%.
That is a price worth paying for a stratum that would otherwise be fitted on noise, and **it is not a
price to pay silently**: every merged stratum carries its `fitted_over`. *Oracle:* two identical
strata merged cost **exactly** nothing, which is the control the harness runs; and a deliberately
non-monotone synthetic sequence must trigger the merge rather than being accepted, while a monotone
one must pass through untouched. *Depends:* E3. *Source:* arch §4.1, spec §4.3, research note §6.6.

**E5. The summary, which is the part a person reads.**  ☐
Several hundred fits per sample against the generic path's four, so the diagnostics **aggregate**
rather than accumulate: how many strata were fitted in place, borrowed and merged, and which; how
many kept their level but borrowed their two shares, which is a different claim from a borrow and is
counted apart from one (spec §4.5); how many fits disagreed across their starting points, with the
worst named; how many carry a guard share above the limit, with the worst named; how many loci stood
behind the thinnest and thickest fit; and `low_slippage_substitution` — the rate this read group's
least-slippery strata returned, which this unit emits and does not compare, because the comparison
needs the generic half.
The per-stratum record is still written — a fit that looks wrong has to be traceable — but nothing
downstream is expected to read it. **A flag nobody reads is how a badly-fitted parameter reaches a
caller**, which is why this step exists at all. *Depends:* E4. *Source:* arch §4.3, spec §4.4.

> **Checkpoint E:** all four parameters are fitted, and each is proven against the harness answer it
> has to reproduce rather than against itself. Pause for review.

### Milestone F — the entry points and end to end

**F1. The two ways in.**  ☐
`SsrEstimationConfig`, `estimate_ssr_parameters(loci, config)` for a caller with nothing else to do
with the stream, and `SsrAccumulators::estimate(config)` for one that drove the accumulator itself.
The first is the second over an accumulator fed by the stream, so the two cannot diverge. A
`LocusGenerationError` is **fatal and propagates**: loci a walk failed to produce are missing
evidence, not zero evidence, and a stratum fitted over a truncated set of loci is wrong in a way
nothing downstream would announce. The read-admission policy travels into the config and out with the
parameters, because a rate describes the reads that survived admission. *Depends:* E5.
*Source:* arch §1.

**F2. Recovery from a directly-filled table — no reads, no reference.**  ☐
Fill a `StratumTable` entry by entry from known parameters and refit, at ploidy 2 **and 4**, and at
three reads a locus **and at forty-five**. The depth arm matters because the scoring rule's exact
unbiasedness was measured to 45 reads a locus and the cap is a precision trade rather than a
correctness limit — a test that only ever runs at three reads cannot tell those apart.

**Fill each entry with its *probability* under the truth rather than with a drawn count**, which is
what makes the tolerance zero: with no sampling noise in the table, "unbiased" is decided rather than
estimated, and anything other than **0.000% on the level and 0.0000 on both shares** is this code's
fault. A table filled by drawing counts instead would need a tolerance nobody could justify, and
would turn the sharpest test in this plan into the vaguest. *Depends:* F1. *Source:* spec §10.1,
§10.2, research note §6.8.

**F3. The identities that need no truth set, on both real cohorts.**  ☐
Three assertions on the tomato CRAMs and the HG002 alignments as they stand, **each run against a
repeat catalog for that reference** (preconditions): one sample walked in one region and in many
gives **identical** tables, which integer entry counts make an equality rather than a tolerance;
`adjustments().loci_without_whole_repeat_reference` is near zero, and a large count is a bug report
against the catalog reader's classification rather than something this unit absorbs; and the table
stays the small object it was measured to be — **tens of thousands of entries, well under one per
locus**, rather than the one-entry-per-locus growth the design would have had to defend.

**The third assertion is a shape and not a number, and that is a change.** The measured walk gave
70,305 entries over 1.73 million tomato loci and 12,727 over 29,811 HG002 loci
([`ng_str_table_memory.rs`](../../../../examples/ng_str_table_memory.rs), research note §6.8) —
under the copy floors of the day, which have since taken nine in ten non-mononucleotide loci off
this path. **Asserting those four numbers would fail on arrival and would say nothing when it did.**
So this step re-runs that program against the real accumulator, asserts entries-per-locus stays far
below one, and **records what it now measures** in the research note, which is where the next size
claim comes from. *Depends:* F2. *Source:* arch §8, research note §6.8, spec §5.1.1.

> **Checkpoint F:** the STR path runs end to end on real alignments, sharded accumulation is an
> equality, and the table's measured size holds against the real implementation. Pause for review.

### Milestone G — the anchor that does not come from the model

**G1. Agreement with HG002's truth genotypes.**  ☐
The parameters fitted by the marginal likelihood must match those measured directly on
known-homozygous loci, **both measured on the same run** rather than one of them read off this page.
**This is the only check in the whole design that does not generate its data from the model it then
fits**: every recovery test above draws from the model, so a shared misspecification cancels and
passes. It is also **the test production's estimator fails**, by 2.4-fold with the direction
reversed, which is the reason for the entire design.

**Measure both sides in this run, and compare them to each other.** The recorded pair — 2.0%
slippage at six or more repeats and a 3.4× direction split (spec §10.3) — was taken before the copy
floors moved, so it is the shape to expect and not the value to assert: a run today fits a shorter
ladder of strata, and a summary pooled over "six or more repeats" is dominated by where its strata
stop (spec §10.3's own trap). **What must hold is that the fit agrees with the direct measurement on
the same loci, and above all that it does not invert the direction split** — reporting gains as more
common than losses is precisely what the biased estimator does. *Depends:* F3. *Source:* spec §10.3,
[`parameter_prepass.md`](../spec/parameter_prepass.md) §2.2.

**G2. The guard share separates strata on real data rather than flagging everything or nothing.**  ☐
Under the old copy floors the walk put 10 of 132 human strata and 28 of 148 tomato strata above the
one-in-ten threshold, so the diagnostic discriminated. **Assert that it still divides the strata —
some above the limit and most below — rather than the two counts**, which the floors have moved:
the strata the floors removed are exactly the low-repeat ones that were failing the threshold, so
the share above it should now be *smaller* and may reach zero on one cohort. **A guard share that
fires on nothing is not conservative, it is absent** — so if it does reach zero, that is a finding
about the floors doing their job and it belongs in the research note, not a green test with nothing
behind it. Also assert that a stratum above the limit is distinguishable in `StratumFitSummary` from
one that merely had few loci. *Depends:* G1. *Source:* spec §5, §10.4, research note §6.8.

**G3. A deliberately unfittable stratum reaches the summary.**  ☐
Feed loci generated with the slippage level at zero and the alleles spread — a stratum whose level is
not identified — and assert it appears in `StratumFitSummary`, named. **A per-stratum record only a
debugger would open does not satisfy this**, which is the whole point of §4.4 and the reason the
summary exists. *Depends:* G2. *Source:* spec §10.6, arch §8.

**G4. A barely-stuttering stratum keeps what it measured and gives up what it did not.**  ☐
Generate a stratum at the 0.091% level §5 measures below four repeats, at a locus count a real
low-repeat stratum has. **Two assertions and they pull in opposite directions:** the level comes back
within its own sampling error — about 5% of itself, which the test computes from the drawn counts
rather than hardcoding — and the direction split and the fall-off come back marked as borrowed, with
`shares_fitted_over` naming their source, rather than as the value about 5 second-step reads happened
to land on. **This needs drawn loci and not the exact method**, unlike everything above it: the whole
question is what sampling does to an estimate whose truth sits against the zero boundary, and the
exact method has no sampling in it. *Depends:* G3. *Source:* spec §4.5, §10.8.

**G5. The two halves of step 4 agree on the substitution rate where nothing slips.**  ☐
**Gated on the generic path's Milestone E1**, and the only step here that is. On one real sample,
this path's substitution rate at its least-slippery strata against the generic path's for the same
read group, on strata at full tract purity: within a quarter-Phred, 6% of the rate. **This is the
second check in the design that does not generate its data from the model it tests** — the other is
G1's — because the two rates are fitted by different models from different sites, so a
misspecification shared with the harness cannot make them agree. **The test lives at step 4's own
surface** (`parameter_estimation/mod.rs`), not in `ssr/`, which never sees the generic half.
**A gap is a result before it is a defect** (spec §1.1): tract impurity charges interruptions to
this path's rate alone, and an error rate genuinely higher inside tracts is the unmeasured finding
§1.1 exists to protect. So the failure message names the purity of the strata compared, and the
limit is not widened to accommodate a gap that survives it. *Depends:* G4. *Source:* spec §4.5,
§8.11, §10.8.

> **Checkpoint G:** the fitted parameters agree with assembly truth on the one dataset that has it,
> the two diagnostics are proven to discriminate rather than to decorate, and a stratum that barely
> stutters is proven to keep its level and disown its shares. Pause for review.

### Milestone H — loci that are not the tract their stratum thinks they are

**Why this is a milestone and not a paragraph in Milestone E** (owner, 2026-08-11). A stratum's four
numbers are pooled over hundreds of thousands of loci on the assumption that each is one tract of
that period and repeat count. Some are not, for biological reasons: a duplication the reference does
not carry collects two paralogous tracts' reads at one locus, a minority of reads mismaps in from a
similar tract elsewhere, a long tract is unstable within the individual, an interruption or a nearby
indel anchors part of the reads differently (spec §4.6). **The sibling path met the same thing and
it inflated the sample's heterozygosity by 41%**: one substitution rate per read group fitted the
body of its distribution and not its tail, and because the only class that could explain such a
locus was the heterozygote, the surplus arrived as heterozygosity. It took a second class of site to
fix, added on its own evidence after that path was already running end to end — which is exactly
where this milestone sits.

**It comes after Milestone G on purpose.** G1 is a real-data comparison against truth genotypes, so
it is the first thing in this plan that could *show* the problem rather than simulate it; and if H
ends in a second class of locus, G1 is one of the checks that re-runs.

**H1. The diagnostic: how much of a stratum sits in shapes its own fit cannot explain.**  ☐
`StratumFit::unexplained_locus_share` and the summary's `worst_unexplained_locus_share` (arch §2.4,
§4.3): the share of a stratum's loci whose shape scores below `UNEXPLAINED_SHAPE_LN_LIMIT` under the
fitted model, **per read and not per locus** — a shape's likelihood is a product over its reads, so
a floor on the total would flag the deepest loci of every stratum by arithmetic. **Free** — every
candidate already evaluates every entry, and an entry is one locus's shape, so this is one more pass
over numbers the search's last step held; no coordinate is kept and the entry key does not move. **It is not the guard share**, which counts *reads* moving by a
non-whole number of copies pooled over a stratum: an aberrant locus can move every read by whole
copies, and one locus in a thousand is invisible in a stratum-wide ratio. Unit tests: a clean
generated stratum reports ~0; the same stratum with one planted locus whose reads sit at two
lengths six repeats apart reports it, and the guard share does not move. Then **read it on tomato
and HG002**, because that is what says which contamination shares H2 has to cover.
*Depends:* G1. *Source:* spec §4.6, arch §2.4.

**H2. What such loci cost, measured exactly.**  ☐ **Own commit, do not bundle.**
In [`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs), generate a stratum
in which a share `q` of the loci come from a **different process** and fit the ordinary model over
it. Two processes, because they fail differently: reads split between two tracts several repeats
apart, which is what a duplication the reference does not carry gives; and reads spread far wider
than any slip kernel, which is instability or mismapping. Report the bias in the level, both shares
and the substitution rate against `q`, over the range H1 measured on real data.

**The exact-bias method works here and that is worth saying, because it does not work for two other
open questions in this plan.** The contaminating process is *known*, so each entry's probability
under the mixture is computable and the table can be weighted by it — no draws, no sampling noise,
so "unbiased" is decided rather than estimated. **The silent failure this isolates:** this path has
two absorbers where the generic path had one. A locus with reads at two distant lengths can be
explained as a locus carrying two distant alleles, which leaves the slippage parameters clean and
**corrupts the emitted allele spectrum instead** — so a run that reports only the four numbers can
pass while the thing it damaged goes unmeasured. *Oracle:* the bias at `q = 0` is exactly zero,
which is the existing control; and the report names **which** parameter moved, the allele
frequencies included. *Depends:* H1. *Source:* spec §4.6, §8.12, §10.9.

**H3. The decision, taken on H2's numbers rather than in advance.**  ☐
If the bias at the shares real data shows sits inside `START_AGREEMENT_LIMIT` — 6%, the finest
difference this design claims a caller can feel — record the measurement, set
`UNEXPLAINED_SHAPE_LN_LIMIT` from where the two populations separated, and stop. If it does not,
**the fix is a second class of locus and nothing else**: ordinary with probability `1 − w`, aberrant
with probability `w` and its reads drawn from a much flatter distribution over the offsets, `w`
fitted rather than chosen. It costs one parameter, it touches nothing shared — the noise parameters
are an associated type of the fitting seam, which is how the generic path added its own second class
without disturbing this one — and it lands the way that one did: **its own milestone, on its own
written-up evidence, with the owner's approval before it is built.**

**What is refused in advance, whatever H2 says: dropping the loci that score badly.** That is
threshold-then-count, the exact bias this whole step exists to remove, and it would take real long
alleles before it took artefacts, since both sit where the model has least mass.
*Depends:* H2. *Source:* spec §4.6, arch §5, §7.

> **Checkpoint H:** what an aberrant locus costs is a number rather than a worry, the parameter that
> absorbs it is named, and the diagnostic that finds such loci on real data has a threshold set by
> measurement. Pause for review — and if a second class of locus is needed, that is a new plan and
> not an extra step here.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | every constrained rate rejects out-of-range and `SsrPeriod` rejects 0 and 7; the allele support's low-end clip pinned by test at repeat counts 3, 6 and 20 |
| B | a table split arbitrarily and merged equals the unsplit one entry for entry in either order; the substitution rate proven to be the maximum by a search that agrees with the closed form to four decimal places (research note §6.7) |
| C | a partial witness proven not to enter as a lost repeat, with the guard removed shown to change the shape; the cap's draw hypergeometric in mean and variance and reproducible from the locus position alone; three shards merged in every order identical, counters included |
| D | the three algebraic gates before any fit — sums to one, no negative counts, silent at a zero level — with the un-rescaled plug-in shown to fail the first at 0.9488; agreement with the harness's kernel to floating point; **the control at exactly 0.000% on the level and 0.0000 on both shares, four starts to 1.000×**, paired with the score at the truth unbeaten |
| E | the harness's own answers: the substitution rate recovered by search, two identical strata merged at exactly zero cost, a non-monotone sequence proven to trigger the merge and a monotone one proven not to |
| F | recovery from a directly-filled table at ploidy 2 and 4 and at 3 and 45 reads a locus, to zero bias; sharded equals single as an equality; the table proven to stay far below one entry per locus on both cohorts, with the current counts recorded |
| G | **the fitted parameters against HG002's known-homozygous measurement, both sides measured in the same run, with the direction split proven not to invert — and, at a low-slippage stratum, this path's substitution rate against the generic path's within a quarter-Phred: the two checks not generated from the model they test**; the guard share proven to divide real strata rather than to fire on all or none; an unfittable stratum proven to reach the summary; a stratum at 0.091% proven to keep its level within sampling error while marking its two shares borrowed |
| H | the bias in all four numbers **and in the emitted allele spectrum** against the share of loci drawn from another process, computed exactly because that process is known, with zero bias at a zero share as the control; the diagnostic proven to see a planted aberrant locus that the guard share cannot |

## Out of scope (next plans)

- **The two censuses** — spec
  [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md); needs an
  architecture document. It is what would make the per-locus stutter model (spec §8.5) and the tied
  genotype-frequency form (spec §8.6) askable at all.
- **The cohort gather** — spec [`parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md);
  needs an architecture document. It owns STR diversity and the reconciliation with this path's
  per-stratum allele spectra.
- **The per-period copy floors** — spec §5.1 and §8.9, blocked on separating the two roles
  `MinCopies` plays in region typing (step 3), not on this step. They change which tracts arrive
  here and nothing about what is done with them.
- **The mismatch count from the aligner** — arch §2.3's `OPEN`, and the better of its two answers.
  Changes one function's body and no signature.
- **A second class of locus, if Milestone H says the fit needs one.** A new plan on H2's written-up
  evidence and with the owner's approval, exactly as the generic path's second class of site was —
  not an extra step bolted onto H3.
- **Three measurable questions this plan does not settle, all cheap and all wanting the harness
  rather than the implementation:** whether a thin stratum should use free genotype frequencies or an
  allele spectrum plus the sample's inbreeding coefficient (spec §8.6); how often the monotonicity
  constraint fires on a truly monotone sequence (spec §8.7); and how a fitted level behaves when the
  truth sits near the zero boundary (spec §8.10). **The last two need random draws rather than the
  exact method, and for the same reason** — a spurious merge and a level piling against zero are both
  sampling effects, and the exact method has no sampling in it. Add worlds to
  [`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs) rather than building
  anything new.
- **What the two substitution rates actually come out at** (spec §8.11). G5 asserts they agree to a
  quarter-Phred, but that limit is borrowed from the generic path's ladder spacing for want of a
  measurement — the first run of the comparison is what sets it, and a persistent gap is §1.1's
  unmeasured finding rather than a failing test.
