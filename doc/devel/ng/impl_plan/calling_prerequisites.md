# ng calling prerequisites — implementation plan

**Status:** draft, 2026-08-21. The build order for **what calling needs from modules it does not
own**: two changes to the cohort merge, three to the parameter pre-pass, and one to `types.rs`.
Each item is already recorded in the settled calling design —
[`read_likelihoods.md`](../spec/read_likelihoods.md) (spec) and
[`../arch/read_likelihoods.md`](../arch/read_likelihoods.md),
[`../arch/calling_priors.md`](../arch/calling_priors.md),
[`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) (types & interfaces) — as a requirement
on an upstream module. This plan turns those recorded requirements into build order; it is **not**
a place for new design. None of it touches frozen production (`src/ssr/`, `src/var_calling/`);
every file changed is ng code.

**Where this sits in the calling build.** Six plans build calling:
`calling_prerequisites` ∥ [`calling_foundations`](calling_foundations.md) →
[`calling_prior`](calling_prior.md) ∥ [`calling_read_likelihoods`](calling_read_likelihoods.md) →
[`calling_loop`](calling_loop.md) → [`calling_bakeoffs`](calling_bakeoffs.md). This plan runs **in
parallel with the foundations plan** and neither needs the other. Downstream the two fan-out plans
need it unevenly: **the read-likelihoods plan needs Milestones B–F here (items 1–5); the prior plan
needs nothing from this plan at all** — its upstream inputs (the fitted spectrum, θ, `F`) already
exist. That asymmetry is why the prior can start the day foundations merges, while the
read-likelihoods plan waits for this one.

---

## Scope

**In:** the six owed items, verbatim from the calling docs:

1. the merge's per-allele support gains a **read-group axis** — one row per `(allele, read group)`;
2. the merge **keeps partial observations** instead of discarding every non-`Complete` witness;
3. the pre-pass emits the **calibration accumulator** — the per-read-group numerator/denominator
   the likelihood's error-rate scale divides;
4. ~~the pre-pass's contamination side-pass emits the contaminating population's three allele-class
   frequencies~~ — **withdrawn 2026-08-24 (owner)**: the mixture's second half is the locus's own
   allele frequency, which the calling loop already estimates, so nothing here is owed. Milestone E
   records why;
5. a **`StratumFits` gather** — the one borrow of `(read group, stratum)` slippage numbers, level
   read off the fitted curve, that crosses the calling seam;
6. **`InbreedingF` tightened to `[0, 1)`**, with the fitted path clamping rather than panicking.

**Out (with owners):**

- **The merge's locus-existence rule for repeat tracts** — counting a partial as non-reference over
  its witnessed stretch so an allele too long to span is not read as "nothing varied"
  ([`read_likelihoods.md`](../spec/read_likelihoods.md) §5.4.2). Owned by whoever brings the STR
  path through the merge; the STR evidence calling consumes today comes straight from the locus
  generator, so nothing in the calling plans is blocked on it.
- **Everything in `src/ng/calling/`** — the four calling plans named above.
- **The choice between the pre-pass's two error-rate routes** — the histogram fit and the census
  fit both exist and which survives is that module's open comparison
  ([`parameter_prepass.md`](../spec/parameter_prepass.md) §4.1). Milestone D gives **both** routes
  an accumulator, exactly because the spec requires the surviving route to carry its own.

## Principles (how the order was chosen)

- **Cheapest-first where order is free.** The six items are independent; the one that shares a file
  with the parallel foundations plan (`types.rs`, Milestone A) lands first so the two branches
  overlap for the shortest time.
- **Types first, then implementation**, within every milestone (project rule).
- **Existing tests are the regression oracle.** Items 1 and 2 change the merge's output shape; a
  sample with one read group must fold to today's shape and today's tests must stay green
  unchanged, which is the parity claim the arch makes ("folding to today's shape where a sample has
  one group", [`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.1).
- **Isolate the silent steps.** A wrong read-group boundary, a mis-projected partial, or an
  accumulator summed over the wrong site set corrupts a downstream genotype without crashing.
  Those steps land as their own commits, oracle green before and after, marked below.
- **Container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); a native host build at
  completion.

## Preconditions (already in place)

- The cohort merge is built ([`cohort_merge.md`](cohort_merge.md) plan complete):
  `SampleSupport`/`SupportedAllele`/`AlleleSupport`
  ([`build.rs:858`](../../../../src/ng/run/cohort_merge/build.rs),
  [`:913`](../../../../src/ng/run/cohort_merge/build.rs),
  [`:973`](../../../../src/ng/run/cohort_merge/build.rs)) and its test suite.
- The pre-pass's two error-rate routes exist: `ReadGroupErrorRateFit`
  ([`generic/read_group_error_rate.rs:45`](../../../../src/ng/parameter_estimation/generic/read_group_error_rate.rs))
  and the joint fit ([`joint/fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs)).
- The contamination estimator exists at the read-group grain: `ContaminationEstimate`
  ([`joint/contamination.rs:430`](../../../../src/ng/parameter_estimation/joint/contamination.rs),
  grain enum [`:238`](../../../../src/ng/parameter_estimation/joint/contamination.rs)).
- The slippage pieces exist: `Slippage`
  ([`joint/ssr_fit.rs:83`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)),
  `StratumFit` ([`:281`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)),
  `blend_level` + `LevelSource`
  ([`joint/slippage_curve.rs:574`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs),
  [`:517`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs)).
- The generic locus generator mints the per-read error the accumulator must match:
  `phred_to_ln_perr(bq).max(mq_log_err)` at
  [`pileup/open_record.rs:2047`](../../../../src/ng/locus_generation/pileup/open_record.rs) and
  [`pileup/fast_column.rs:211`](../../../../src/ng/locus_generation/pileup/fast_column.rs).

## Worktree, branch, merge

- **Worktree** `../pop_var_caller-calling-prerequisites`, **branch** `ng-calling-prerequisites`,
  from `main`, plain `git worktree add` (repo convention).
- **Runs in parallel with** `ng-calling-foundations`. The shared file is `src/ng/types.rs`, and the
  overlap is avoided by region, not resolved by merge: this branch edits **only** the `InbreedingF`
  block ([`types.rs:388`](../../../../src/ng/types.rs)), its boundary test
  ([`types.rs:862`](../../../../src/ng/types.rs)), and inserts its new `DomainError` variant
  **immediately after the existing `InbreedingF` variant**; foundations appends its new scalars and
  variants **at the end** of their sections. Disjoint regions, no textual conflict expected.
- **Merge order back:** whichever of the two finishes first merges to `main` first; the second
  merges `main` in and re-runs its tests. If the `DomainError` enum does conflict anyway, the
  second-merger resolves it (it is a variant append on both sides).
- Milestones B–F touch `cohort_merge/` and `parameter_estimation/`, which no other calling branch
  edits.

---

## The steps

### Milestone A — `InbreedingF` in `[0, 1)`, with the clamp (item 6)

The ceiling is a property of the type, not of one estimator
([`calling_priors.md`](../spec/calling_priors.md) §7). The arch names the three-part blast radius
([`../arch/calling_priors.md`](../arch/calling_priors.md) §2.1), and each part is a step.

**A1. The half-open check.**  ✅
`InbreedingF::try_new` rejects `1.0`: its own `[0, 1)` range test with a new `DomainError` variant
that says so, **not** a change to the shared `checked_probability`
([`types.rs:326`](../../../../src/ng/types.rs)), which the other fraction newtypes share and which
is right to admit `1.0` for them. Move the existing acceptance assertion at
[`types.rs:862`](../../../../src/ng/types.rs) to the rejection list beside `1.5`. *Source:*
calling_priors arch §2.1; spec §7.

**A2. The fitted path clamps instead of panicking.**  ✅
[`runs.rs:634`](../../../../src/ng/parameter_estimation/generic/runs.rs) builds an `InbreedingF`
from a coverage-weighted posterior occupancy with `.expect(…)`; that occupancy can in principle
reach exactly `1.0` on a fully homozygous sample, so after A1 the `expect` is a panic on a
legitimate fit. Replace it: clamp the fitted value at **`0.99`** before constructing — production's
own estimator ceiling, imported with its reasoning ("no sample ever reaches the caller carrying a
prior that has ruled heterozygotes out",
[`paralog/inbreeding.rs:25`](../../../../src/paralog/inbreeding.rs); spec §7). Test: a fit of
exactly `1.0` constructs at `0.99` and does not panic. **Own commit, do not bundle** — a wrong
clamp is a quietly different prior, not a crash; the oracle is A1's boundary tests plus the
clamp test. *Depends:* A1. *Source:* calling_priors arch §2.1; spec §7.

> **Checkpoint A:** `InbreedingF` rejects `1.0`, every constructor site compiles, the fitted path
> clamps, `cargo test` green. Pause for review.

### Milestone B — the merge's read-group axis (item 1)

Summing must stop at the read-group boundary: two reads showing the same sequence from two lanes
have different error rates and must not be pooled ([`read_likelihoods.md`](../spec/read_likelihoods.md)
§2.3). The merge's own doc already books the change as owed
([`build.rs:958`](../../../../src/ng/run/cohort_merge/build.rs)).

**B1. The type change.**  ✅
`SupportedAllele` ([`build.rs:913`](../../../../src/ng/run/cohort_merge/build.rs)) gains the
read group: one row per `(allele, read group)`, rows in ascending `(allele, read group)` order —
the shape that folds to today's where a sample has one group. `ReadGroupId` is already on
`SequenceObservation` ([`locus_generation/mod.rs:316`](../../../../src/ng/locus_generation/mod.rs)).
Compile-driven follow-through on every consumer of `supported`. *Source:* read_likelihoods spec
§2.3; arch §2.1.

**B2. Attribution stops at the boundary.**  ✅
The collation that today merges a sample's observations "where two of its own observations reached
the same allele" now merges only within one read group; the divided-read paths
(`AlleleSupportTally`) key their tallies by `(allele, read group)`. Tests: a two-read-group fixture
keeps two rows for one allele, with `num_reads` and `q_sum` split exactly as the reads were; a
one-group fixture is **byte-identical to today's output** — the existing merge test suite is the
oracle and must pass unchanged apart from the row type. **Own commit, do not bundle** — pooling
across the boundary is a quietly-wrong `q_sum`, not a crash. *Depends:* B1. *Source:*
read_likelihoods spec §2.3; arch §2.1; [`build.rs:958`](../../../../src/ng/run/cohort_merge/build.rs).

> **Checkpoint B:** the axis exists, one-group samples are unchanged, two-group samples split.
> Pause for review.

### Milestone C — the merge keeps partial observations (item 2)

Today collation skips every observation whose witness is not `Complete`
([`build.rs:1351`](../../../../src/ng/run/cohort_merge/build.rs)) and projection panics rather than
pad one ([`build.rs:323`](../../../../src/ng/run/cohort_merge/build.rs)), so the likelihood's
censored term has nothing to read. The requirement, verbatim from the corrected spec: **a partial
observation must survive collation, keyed by the stretch it witnessed, and projected over that
stretch rather than the whole locus span**
([`read_likelihoods.md`](../spec/read_likelihoods.md) §5.4, corrected 2026-08-21; §7's ownership
table). It changes **no locus's existence** — the variability filter still counts complete
observations only on this path (§5.4.2).

**C1. The carried type.**  ✅
A partial row on `SampleSupport`: the witnessed stretch (offset + length, off `ReadWitness`), the
bases over that stretch, `num_reads`, `q_sum`, and the read group (B's axis applies here too —
the evidence view's `PartialObservation` is consumed per read group like everything else).
No logic. *Depends:* B1. *Source:* read_likelihoods spec §5.4; arch §2.1 (`SampleEvidence.partials`).

**C2. Collation keeps them; projection projects the stretch.**  ✅
The `!= Complete → continue` at [`build.rs:1351`](../../../../src/ng/run/cohort_merge/build.rs)
routes partials into C1's rows instead of dropping them; the projection gains a
witnessed-stretch variant instead of the panic (the panic stays for the code path that must never
see one). Tests: a fixture with one partial read yields the row with hand-written stretch and
bases; a fixture with none is byte-identical to today (existing suite green); locus existence
verdicts unchanged on every existing fixture. **Own commit, do not bundle** — a partial projected
over the whole span mis-scores as a short allele silently
([`read_likelihoods.md`](../spec/read_likelihoods.md) §5.1); the oracle is the hand-written
fixture plus the untouched existing suite. *Depends:* C1. *Source:* read_likelihoods spec §5.4,
§5.1.

> **Checkpoint C:** partials survive the merge, keyed and projected over their stretch; nothing
> else about the merge's output moved. Pause for review.

### Milestone D — the calibration accumulator (item 3)

The likelihood's scale is `fitted rate / mean minted per-read error`
([`read_likelihoods.md`](../spec/read_likelihoods.md) §3.2). Two scalars per read group — the sum
of minted per-read error probabilities and the count of reads summed — with two requirements: the
minted quantity is **the same function the locus generator mints with**, and the sum runs over
**exactly the sites the surviving error-rate estimate was fitted from**, per route.

**D1. One mint function.**  ✅
Hoist the per-read error mint — worse of the window's base quality and the mapping quality, in log
space — into one named `pub(crate)` function beside its current home, and call it from both mint
sites ([`pileup/open_record.rs:2047`](../../../../src/ng/locus_generation/pileup/open_record.rs),
[`pileup/fast_column.rs:211`](../../../../src/ng/locus_generation/pileup/fast_column.rs)). Test:
byte-identical `q_sum` on the existing pileup fixtures. *Source:* read_likelihoods spec §3.2,
§12 test 10 ("checked by calling it from both sides on the same read").

**D2. The accumulator, on the route that can carry it.**  ✅
**The step as written asked for something that does not exist, and the owner settled it on
2026-08-24. Two corrections, both now in
[`read_likelihoods.md`](../spec/read_likelihoods.md) §3.2.**

**The average is the geometric mean, not the arithmetic one.** §3.2 asked for a running sum of the
per-read error *probability*, and nothing carries that: the walk sums the *logarithms* into an
observation's `q_sum` and throws the individual reads away, so `Σ ε` is not recoverable. Supplying it
would have meant a second accumulation at fold time and a new field on every observation. Taking
`exp(Σ q_sum / Σ num_obs)` instead is the self-consistent choice rather than a concession: it is the
quantity the model charges an observation, and the one production uses
([`posterior_engine.rs:1536`](../../../../src/var_calling/posterior_engine.rs) — production has no
recalibration at all, so there is nothing there to copy but the quantity). **So the accumulator adds
up numbers that already exist**, which is what "no new traversal" meant all along.

**Neither route can call D1's function, and neither needs to.** The census route's per-position unit
is a depth code and allele counts ([`joint/fit.rs:467`](../../../../src/ng/parameter_estimation/joint/fit.rs))
with no quality in it, and the histogram route reads pooled observations rather than reads. The
histogram route can supply both numbers from the observations it already walks; **the census route
cannot supply either, whichever average is chosen**, so its accumulator waits on §4.1's comparison
between the two routes — if it wins, its records gain a quality field; if the histogram route wins,
nothing is owed.

**Built:** [`generic/calibration.rs`](../../../../src/ng/parameter_estimation/generic/calibration.rs),
accumulated in `GenericAccumulators::add_locus` over exactly the loci and read groups the
error-rate histogram counts, exposed beside the histograms because nothing is fitted from it.
**The sum is a fixed-point integer**, because the accumulator promises order-independent merging
across region shards and an `f64` running sum would have broken that silently. *Source:*
read_likelihoods spec §3.2; arch §3.

> **Checkpoint D:** the histogram route carries the accumulator, over the sites its own fit reads,
> minted by nothing — the numbers come from the walk. The census route's is deferred to the
> comparison that decides whether that route survives at all. Pause for review.
>
> **Reviewed 2026-08-24**, four ways, in
> [`ng_prereq_closeout_d2_review`](../../reports/implementations/ng_prereq_closeout_d2_review_2026-08-24.md)
> — the step had had none. Six defects, all fixed here except one, which is the owner's: **the
> per-position depth cap does divide the two site sets**, by 2.7% on HG002 at 300× and by nothing on
> tomato, where the spec said it does not. Both options are in
> [`read_likelihoods.md`](../spec/read_likelihoods.md) §3.2. And **the two averages the choice of
> mean is between are 25 to 44 times apart** on real reads, not close as the decision assumed —
> [`ng_prereq_closeout_two_averages`](../../reports/implementations/ng_prereq_closeout_two_averages_2026-08-24.md).

### Milestone E — deleted: the contamination mixture's second half is not the pre-pass's (item 4)

**This milestone asked the parameter pre-pass for the contaminating population's three allele-class
frequencies — how often it carries the reference, a substitution, an insertion or deletion — from a
side-pass over the census sites. That is deleted rather than built (owner, 2026-08-24), and nothing
replaces it here.**

**What the model needs is the frequency of the allele an observation shows, at the locus being
called** ([`read_likelihoods.md`](../spec/read_likelihoods.md) §3.6, corrected the same day). The
three-class split is not a modelling requirement: it is what production does because it has no
per-locus population frequency for an arbitrary alternative allele and must fall back on a class
average. This caller has that frequency — it is the same one the genotype prior reads and the
calling loop re-estimates — so the classes disappear along with the ignorance that motivated them.

**Three reasons this could never have been the pre-pass's**, all of which came out of trying:

- **The pre-pass visits a selection of census sites; the caller calls everywhere.** A side-pass
  cannot supply a per-locus frequency for a locus it never looked at, so no amount of work here
  would have produced the number.
- **The estimator could not have been ported.** Production fits its three-entry simplex from
  per-read posteriors that a read came from the contaminant
  ([`var_calling/contamination_estimation.rs`](../../../../src/var_calling/contamination_estimation.rs)).
  ng's contamination fit computes no such posterior — it works from marker read shares against
  ancestry-predicted dosages — so there was nothing to take shares of, and inventing an estimator is
  design this plan says in its own header it is not a place for.
- **The census cannot answer the indel half cleanly anyway.** Its fifth allele code lumps insertions
  and deletions together with `N` and spanning deletions
  ([`joint/census.rs`](../../../../src/ng/parameter_estimation/joint/census.rs),
  `ObservedAllele::Other`), so an "insertion-or-deletion frequency" read off it is contaminated by
  two things that are not alleles.

**Owner:** the read likelihood, in
[`calling_read_likelihoods.md`](calling_read_likelihoods.md). It is a few lines there because the
frequency is already in hand. The one consequence recorded with the design: the contamination
*fraction* stays frozen before the loop and this frequency does not, so the two halves of the
mixture sit in different tiers.

### Milestone F — the `StratumFits` gather (item 5)

The calling seam needs one borrow: the `(read group, stratum)` slippage lookup with the **level
read off the fitted curve** rather than the cell
([`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §4.2). The pieces exist; nothing
gathers them, and the loop's `FrozenParameters` sketch names the missing wrapper `StratumFits`
([`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §2, open items).

**F1. The gather.**  ✅
`StratumFits` in `parameter_estimation/joint/`: built once per run from the `StratumFit`s and the
fitted curves; lookup by `(read group, period, repeat count)` returns the stratum's `Slippage`
with `level` replaced by `blend_level`'s value for that cell and the `LevelSource` provenance
carried. Test: the returned level equals calling `blend_level` directly (the curve is the oracle),
the shape numbers equal the stratum's own, and provenance survives. Update
[`../arch/parameter_prepass_joint_fit.md`](../arch/parameter_prepass_joint_fit.md) to name the
type, per the em-loop arch's instruction that the parameter-prepass arch doc pins it. *Source:*
read_likelihoods arch §4.2; calling_em_loop arch §2.

> **Checkpoint E/F:** the gather exists with its tests and the pre-pass's own outputs are
> otherwise untouched. **Every precondition the read-likelihoods plan needs from here now holds or
> has been withdrawn:** the merge's read-group axis, its partial observations and the slippage
> gather are built; the calibration accumulator is built on the histogram route, with the census
> route's deferred to the comparison that decides whether that route survives; and the contamination
> mixture's second half turned out not to be the pre-pass's at all. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | boundary tests both directions; the fitted-`1.0` clamp test |
| B | **parity:** one-read-group fixtures byte-identical to today (existing merge suite); two-group fixture splits rows exactly |
| C | hand-written partial fixture (stretch, bases, counts); existing suite green; locus-existence verdicts unchanged |
| D | same-function check from both mint sites; hand-computed fixture mean; merge in any shard order gives one answer; **on real reads, the totals count exactly the reads the walk emitted at the loci the histogram counts** — 172,616,054 on HG002 at 300×, both paths, `examples/ng_minted_error_means.rs`. *(Not "equals the histogram's own read count": the histogram thins every position to `MAX_BINNED_DEPTH` = 124 and these totals thin nothing, so above 124 reads a position the two numbers differ by design — §3.2's argument is that the cap is quality-blind, so the **mean** is unbiased, not that the counts match.)* |
| E | deleted — nothing to verify here |
| F | **the curve as oracle:** gathered level ≡ `blend_level` called directly; provenance carried |

## Out of scope (next plans)

- **The evidence views over the changed merge types** (`GenericSampleEvidence`,
  `PartialObservation`) — [`calling_read_likelihoods.md`](calling_read_likelihoods.md), which
  consumes what B and C produce.
- **`FrozenParameters` assembling the gather, the calibration and the contamination views** —
  [`calling_loop.md`](calling_loop.md).
- **The contamination mixture's second half** — the frequency of the observed allele among the
  samples in this one's sequencing batch, recomputed each iteration
  ([`calling_read_likelihoods.md`](calling_read_likelihoods.md); spec §3.6, and Milestone E for why
  it left this plan).
- **The repeat-path locus-existence amendment** — whoever brings the STR path through the merge
  (see Scope).
- **Choosing between the two error-rate routes** — the pre-pass's own comparison. D gives the
  histogram route its accumulator; the census route cannot carry one without a quality field in its
  records, so that comparison now also decides whether those records change.
