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
4. the pre-pass's contamination side-pass emits the **contaminating population's three allele-class
   frequencies**;
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

**C2. Collation keeps them; projection projects the stretch.**  ☐
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

**D1. One mint function.**  ☐
Hoist the per-read error mint — worse of the window's base quality and the mapping quality, in log
space — into one named `pub(crate)` function beside its current home, and call it from both mint
sites ([`pileup/open_record.rs:2047`](../../../../src/ng/locus_generation/pileup/open_record.rs),
[`pileup/fast_column.rs:211`](../../../../src/ng/locus_generation/pileup/fast_column.rs)). Test:
byte-identical `q_sum` on the existing pileup fixtures. *Source:* read_likelihoods spec §3.2,
§12 test 10 ("checked by calling it from both sides on the same read").

**D2. The accumulator, both routes.**  ☐
`ReadGroupErrorRateFit`
([`read_group_error_rate.rs:45`](../../../../src/ng/parameter_estimation/generic/read_group_error_rate.rs))
gains the two fields, summed over the sites its histogram counts; the census route's fit
([`joint/fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs)) gains its own pair,
summed over the census sites its fit reads, calling D1's function. Whichever route the pre-pass's
§4.1 comparison keeps, its accumulator goes with it and the other is deleted with the other fit.
Test: on a fixture, `Σ minted / count` equals a hand-computed mean; a route's accumulator counts
no site the route's fit did not read. **Own commit, do not bundle** — numerator and denominator
from different site sets "is not a calibration" and nothing crashes; the oracle is the
hand-computed fixture mean. *Depends:* D1. *Source:* read_likelihoods spec §3.2; arch §3.

> **Checkpoint D:** both routes carry the accumulator, minted by the one shared function. Pause
> for review.

### Milestone E — the contaminating population's allele-class frequencies (item 4)

The mixture's second half: how often the contaminating population carries the reference, a
substitution alternative, an insertion-or-deletion alternative — one frequency per allele class,
averaged over the census sites, following production and freebayes
([`read_likelihoods.md`](../spec/read_likelihoods.md) §3.6).

**Only this half is owed.** The contamination *fraction* and its evidence counts landed on `main`
with the read-group-grain work: `ContaminationEstimate::Estimated` already carries `alpha`,
`source`, `markers_with_reads` and `reads_on_markers`
([`joint/contamination.rs:430`](../../../../src/ng/parameter_estimation/joint/contamination.rs)),
which is four of `ContaminationView`'s six fields. The three class frequencies are the two that
are missing, and `joint/contamination.rs` has none of them.

**The grain question is settled by that same work, not left open.** The spec's §6 table puts the
fraction and the population's frequencies on one row at read-group grain; only the first belongs
there, and the merged code says why in its own words — the fraction is per read group because
*"a second plant's DNA enters at library preparation or at sequencing, so two libraries of one
plant can carry different amounts of it"*
([`:238`](../../../../src/ng/parameter_estimation/joint/contamination.rs)). That reasoning is
about **how much** contaminant there is. The class frequencies say **what** the contaminant
carries, and a run's contaminant is one population however many libraries it entered through. So
they are fitted once per run and copied into each read group's view.

**E1. The three frequencies.**  ☐
A side-pass over the census evidence the contamination fit already reads: classify each site's
alternative alleles into the three classes, average the fitted frequencies per class, and carry
the triple beside the run's contamination output (one triple per run — the frequencies describe
the contaminating population, not any one read group; the per-read-group `ContaminationView`
copies them in, [`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.3). Emitted as
absent where contamination itself is not estimated (one sample). Test: a fixture census with
hand-classified sites gives the hand-computed averages. *Source:* read_likelihoods spec §3.6;
arch §2.3.

### Milestone F — the `StratumFits` gather (item 5)

The calling seam needs one borrow: the `(read group, stratum)` slippage lookup with the **level
read off the fitted curve** rather than the cell
([`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §4.2). The pieces exist; nothing
gathers them, and the loop's `FrozenParameters` sketch names the missing wrapper `StratumFits`
([`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §2, open items).

**F1. The gather.**  ☐
`StratumFits` in `parameter_estimation/joint/`: built once per run from the `StratumFit`s and the
fitted curves; lookup by `(read group, period, repeat count)` returns the stratum's `Slippage`
with `level` replaced by `blend_level`'s value for that cell and the `LevelSource` provenance
carried. Test: the returned level equals calling `blend_level` directly (the curve is the oracle),
the shape numbers equal the stratum's own, and provenance survives. Update
[`../arch/parameter_prepass_joint_fit.md`](../arch/parameter_prepass_joint_fit.md) to name the
type, per the em-loop arch's instruction that the parameter-prepass arch doc pins it. *Source:*
read_likelihoods arch §4.2; calling_em_loop arch §2.

> **Checkpoint E/F:** the class frequencies and the gather exist with their tests; the pre-pass's
> own outputs are otherwise untouched. **This plan is complete — the read-likelihoods plan's
> preconditions (items 1–5) all hold.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | boundary tests both directions; the fitted-`1.0` clamp test |
| B | **parity:** one-read-group fixtures byte-identical to today (existing merge suite); two-group fixture splits rows exactly |
| C | hand-written partial fixture (stretch, bases, counts); existing suite green; locus-existence verdicts unchanged |
| D | same-function check from both mint sites; hand-computed fixture mean; per-route site-set discipline |
| E | hand-classified fixture census → hand-computed class averages |
| F | **the curve as oracle:** gathered level ≡ `blend_level` called directly; provenance carried |

## Out of scope (next plans)

- **The evidence views over the changed merge types** (`GenericSampleEvidence`,
  `PartialObservation`) — [`calling_read_likelihoods.md`](calling_read_likelihoods.md), which
  consumes what B and C produce.
- **`FrozenParameters` assembling the gather, the calibration and the contamination views** —
  [`calling_loop.md`](calling_loop.md).
- **The repeat-path locus-existence amendment** — whoever brings the STR path through the merge
  (see Scope).
- **Choosing between the two error-rate routes** — the pre-pass's own comparison; D leaves both
  routes carrying their accumulator so the choice loses nothing.
