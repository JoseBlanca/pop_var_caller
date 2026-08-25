# ng calling loop (step 9, arm A) — implementation plan

**Status:** draft, 2026-08-21. The build order for **the default calling loop**: the shared
per-locus types the fan-out plans deferred (`CallingScratch`, `LocusEvidence`,
`FrozenParameters`), the `LocusGenotyper` seam, and `SummariseConditionLoop` — arm A, with the
outer two loops **structurally present and switched off**, at the default configuration. Design
is settled in [`calling_em_loop.md`](../spec/calling_em_loop.md) (spec) and
[`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) (types & interfaces). This plan turns
that design into build order; it is **not** a place for new design.

**This plan ends with ng calling genotypes.** Everything after it is measurement: arms B/C, the
exhaustive scorer, the slippage re-fit and discovery mechanisms, and every bench run live in
[`calling_bakeoffs.md`](calling_bakeoffs.md).

**Where this sits.** Six plans build calling:
[`calling_prerequisites`](calling_prerequisites.md) ∥
[`calling_foundations`](calling_foundations.md) → [`calling_prior`](calling_prior.md) ∥
[`calling_read_likelihoods`](calling_read_likelihoods.md) → `calling_loop` →
[`calling_bakeoffs`](calling_bakeoffs.md). This plan needs **both** fan-out plans: it drives the
prior's builders and the likelihood's rows, and owns nothing of their mathematics.

**That blocker is lifted on the SNP/indel path, and this paragraph replaces the one that said
otherwise** (2026-08-25). This plan was written when candidate selection had no spec, and recorded
as its one blocker that every integration test would take a fixture-supplied candidate set and
that an end-to-end run on real data was impossible. **Step 6 now has a spec, an architecture and a
shipped implementation on the generic path**: `select_generic`
([`../spec/candidate_alleles.md`](../spec/candidate_alleles.md),
[`../arch/candidate_alleles.md`](../arch/candidate_alleles.md),
[`candidate_alleles.md`](candidate_alleles.md)) takes one assembled `CohortObservation` and
returns the narrowed `CandidateAlleles` with the per-sample leftover the read likelihood needs.

**What that changes here, and what it does not.** Hand-built candidate sets remain right for the
unit tests, and for a reason stronger than convenience: a test that calls selection and the loop
together cannot say which of the two is wrong. What changes is Milestone E — an integration test on
the generic path may now build its candidates by calling selection, and the end-to-end run on real
data is no longer blocked for that path.

**It is still blocked for repeat tracts, by two independent gaps.** The STR read-likelihood row
does not exist — `censored_emission` is `unimplemented!()`
([`calling_read_likelihoods.md`](calling_read_likelihoods.md) G1; the row is its H1 and H2) — and
**that one blocks a tract genotype outright**, whatever candidates are supplied. Separately, the
STR candidate path is unwritten: `allele_candidates/` holds `mod.rs` and `generic.rs` only, and
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) is a later session's. The first must close
before a tract can be called at all; the second only decides whether its candidates are chosen or
handed in.

---

## Scope

**In:** `calling/mod.rs` gains `CallingScratch`, `LocusEvidence`, `FrozenParameters`;
`calling/inference/mod.rs` — `LocusGenotyper`, `CallingLoopConfig` (with `SlippageRefitConfig`
and `DiscoveryConfig` as **values**, defaults `rounds = 0` / `Off`);
`calling/inference/summarise_condition.rs` — the three nested loops with the outer two inert at
the defaults; the input edge (evidence views built from the merge's output and the locus
generator's; `FrozenParameters` assembled from the pre-pass); the two loop oracles.

**Out (later plans):**

- **Arms B and C, the exhaustive scorer, the local search, `assignment.rs`** —
  [`calling_bakeoffs.md`](calling_bakeoffs.md) (spec Q1).
- **The slippage re-fit body and `discovery.rs`** — the same plan (spec Q2, Q3). Here their
  configs exist and their non-default values are **rejected with a clear message** until the
  bodies land, so a configuration cannot silently no-op.
- **Candidate selection (step 6)** — **built and merged on the generic path; consume it, do not
  edit it.** `src/ng/calling/allele_candidates/` is another plan's
  ([`candidate_alleles.md`](candidate_alleles.md)) and is still moving through that plan's
  Milestone D. Its public surface is settled; the one value that may still change is
  `DEFAULT_MAX_CANDIDATE_ALLELES`. The repeat-tract half is not written at all.
- **Emission, QUAL, site filters, phasing** — steps 10–13
  ([`ng_proposal.md`](../spec/ng_proposal.md)).
- **Where the loop runs inside the merge's builder** — the wiring into `run/` follows the
  end-to-end blocker; this plan's driver takes a `CohortObservation` and is callable from the
  builder when selection exists (spec §9 says the placement commutes).

## Principles (how the order was chosen)

- **Types first, then implementation** (project rule) — the seam and scratch before any pass.
- **The algorithmic heart before the plumbing.** One E/M pass is built and proven on hand-built
  rows before the loop repeats it; the loop is proven before any real evidence is shaped into it.
- **Reuse over rewrite.** The loop's shape, the flat first pass, the convergence arithmetic, the
  emit-with-flag rule and the scratch lift are ports from
  [`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs) (spec §10's map); ng
  needs no inner trait — the paths differ only in the sibling row builders.
- **Verify against ground truth.** The SNP/indel side has a **parity oracle** — production's
  loop on the same likelihood table, prior parameters and candidates must yield the same
  genotypes, any difference traced to a recorded decision. The STR side gets a **differential,
  not a parity oracle**: production's STR loop converges on π at `1e-6`
  ([`em.rs:137`](../../../../src/ssr/cohort/em.rs)) where ng converges on copies over
  chromosomes at `1e-3`, so ng runs once under production's rule (genotypes must match), then
  under its own, and reports what moved — a check with a failing state (spec §10).
- **Isolate the silent steps.** The fixed-order M-step sum, the ÷chromosomes in the convergence
  delta, and the emission-call count are the named quietly-wrong candidates; each is its own
  commit, marked below.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place)

- **The prior plan merged:** `sample_concentration`, `GenotypePriorModel` with both impls, the
  seed builders ([`calling_prior.md`](calling_prior.md)).
- **The read-likelihoods plan merged:** both row functions, the emission seam, the evidence
  views, `ReadGroupCalibration`/`ContaminationView`
  ([`calling_read_likelihoods.md`](calling_read_likelihoods.md)).
- Through those, prerequisites and foundations: the genotype table, the merge's axis + partials,
  the calibration accumulator, the contamination class frequencies, the `StratumFits` gather.
- Production's loop internals as reuse targets: `GenotypeEmModel`/`run_em_loop`/`EmStepPhase`
  ([`posterior_engine.rs:2635`](../../../../src/var_calling/posterior_engine.rs),
  [`:2733`](../../../../src/var_calling/posterior_engine.rs),
  [`:2586`](../../../../src/var_calling/posterior_engine.rs)), the convergence arithmetic and
  comment ([`:2704–2726`](../../../../src/var_calling/posterior_engine.rs)), the constants
  ([`:86`](../../../../src/var_calling/posterior_engine.rs),
  [`:96`](../../../../src/var_calling/posterior_engine.rs)), the retired-error rule
  ([`:26`](../../../../src/var_calling/posterior_engine.rs)), the STR final pass
  ([`em.rs:857`](../../../../src/ssr/cohort/em.rs)).
- **Candidate selection merged**, on the generic path: `select_generic` and its config, verdict,
  remapping and per-sample leftover ([`candidate_alleles.md`](candidate_alleles.md)).
- **Not in place:** the two repeat-tract halves — the STR read-likelihood row and the STR
  candidate path. Both are named in the blocker note above.

## Branch and merge (sequential — no worktree)

- **Branch** `ng-calling-loop`, from `main` **after both `ng-calling-prior` and
  `ng-calling-read-likelihoods` have merged**. Nothing runs beside this plan, so it gets **no
  worktree** — a plain branch in the primary checkout is the convention for sequential work.
- Conflict surface: none expected — every phase-2 branch is already in `main`, and this branch
  is the sole editor of `calling/mod.rs` and `calling/inference/` while it lives.
- Merges straight back to `main`; [`calling_bakeoffs.md`](calling_bakeoffs.md) branches after it.

---

## The steps

### Milestone A — the deferred shared types + the seam (types, no logic)

**A1. `CallingScratch`, `LocusEvidence`, `FrozenParameters`.**  ✅
In `calling/mod.rs`, now that every borrowed field exists: `CallingScratch` (the `Lg` table,
posterior row, concentration, current/previous expected copies, per-sample copies, and the
likelihood's `RowScratch` section — allocated once per worker, reused per locus; the measured
16%-of-cycles allocator reason travels in the doc comment); `LocusEvidence`
(`Generic`/`Ssr` per-sample views — the one place the paths meet); `FrozenParameters`
(calibration, contamination, per-sample `InbreedingF`, `SpectrumSeed`, `&StratumFits`, ploidy).

**`CallingScratch` also gains candidate selection's buffers as a field** — `SelectionScratch`
(`src/ng/calling/allele_candidates/mod.rs`), which that module built standing alone and whose own
doc comment names this step as where it moves. The same worker runs selection and then the loop on
the same locus, so a second per-worker allocation would buy nothing
([`../arch/candidate_alleles.md`](../arch/candidate_alleles.md) §2.4). **Nothing about its shape
changes** — move it in, do not restate it.
The ordering contract documented and asserted: one run-wide sample order indexes every per-sample
slice; `LocusEvidence`'s discriminant and `CandidateAlleles.kind` must agree — disagreement is a
caller bug and asserts. *Source:* arch §2; spec §8.

**A2. `LocusGenotyper` + `CallingLoopConfig`.**  ☐
`inference/mod.rs`: the seam trait (`call_locus(evidence, parameters, candidates, config,
scratch) -> LocusInference`); `CallingLoopConfig` with the inherited constants as named,
soft-marked values — `DEFAULT_CONVERGENCE_THRESHOLD = 1e-3`, `DEFAULT_MAX_PASSES = 50` — plus
`SlippageRefitConfig` (default `max_rounds = 0`) and
`DiscoveryConfig` (default `Off`) as **values, not code paths**. Until
[`calling_bakeoffs.md`](calling_bakeoffs.md) lands the bodies, a non-default setting of either is
rejected loudly at config validation — never silently ignored.

**The allele cap is deliberately not among those constants, and this entry used to declare it**
(corrected 2026-08-25). `DEFAULT_MAX_CANDIDATE_ALLELES` already exists, in
`src/ng/calling/allele_candidates/mod.rs`, as a `MaxCandidateAlleles` — a newtype that refuses
anything below two, because a cap of 0 or 1 is refusal under another name. A second constant of
the same name here, typed `u16`, would be two spellings of one rule and would drop the check;
selection's `CandidateSelectionConfig` is what carries the cap and the support bar, and the loop
takes it rather than restating it. *Depends:* A1. *Source:* arch
§2.1, §3.1, §6.1, §6.2.

> **Checkpoint A:** the seam compiles against both siblings' types; the ordering contract is
> asserted. Pause for review.

### Milestone B — one E/M pass (the heart, on hand-built rows)

**B1. The E-step for one sample.**  ☐
Given a filled `Lg` row: build the sample's concentration (`sample_concentration`), its log-prior
row (`GenotypePriorModel`), add, softmax into the posterior row, fold into per-sample expected
copies. Pure, scratch-backed, no allocation. Test: a hand-computed 2-allele diploid case, prior
and likelihood chosen so every intermediate is checkable by hand. *Depends:* A2. *Source:* spec
§2; arch §1.

**B2. The M-step, fixed order.**  ☐
The cohort's expected copies as a sum over samples **in the run's fixed sample order**. Tests
(spec §13 test 2): permuting the samples changes no genotype; **the mutation check is on the
summed expected copies, compared bitwise** — not on the argmax, which will not flip on a
last-bit reorder and proves nothing. **Own commit, do not bundle** — a reordered float sum is
quietly different output at another worker count, never a crash; the bitwise test is the oracle.
*Depends:* B1. *Source:* spec §8; arch §4.

> **Checkpoint B:** one pass is exact, deterministic, and hand-verified. Pause for review.

### Milestone C — the frequency loop

**C1. The flat first pass.**  ☐
The prior-free initialisation (reads only), ported with its reasoning — `EmStepPhase`'s shape
([`posterior_engine.rs:2586`](../../../../src/var_calling/posterior_engine.rs)); it runs at the
start of every outer round, which the skeleton in D preserves. Test (spec §13 test 3): after the
first pass the expected copies reflect the reads and not the seed; **and the trap itself** — at
a locus thin enough that the seed outweighs each sample's reads, a seeded first pass converges
to no-variant where the flat start converges to the variant. *Depends:* B2. *Source:* spec §3.

**C2. Convergence, the cap, and the emitted flag.**  ☐
Stop when the largest change in expected copies **divided by cohort chromosomes** falls below the
threshold; at the cap, emit with `converged = false` — never dropped, never an error — and
`passes` always emitted. Tests: spec §13 test 4 (capped locus emitted, flag reaches the output;
**no** assertion that the delta falls every pass — EM guarantees a monotone likelihood, not a
monotone delta); the one-sample fixed point (test 1 — pass 2's copies equal pass 1's **bit for
bit**, the loop stops, asserted on the genotype and on pass-1-equals-pass-2, with **no branch on
cohort size anywhere**); and the division's own test — the same raw-copies movement at `n = 1`
and `n = 1000` crosses the threshold at the same *frequency-scale* point, which fails if anyone
drops the divisor. **Own commit, do not bundle** — the ÷chromosomes is load-bearing across the
cohort range and a criterion written on raw counts tightens silently with cohort size; the
two-cohort-size test is the oracle. *Depends:* C1. *Source:* spec §6, §7; arch §4.

**C3. The final pass.**  ☐
Score every sample once more, take the highest-posterior genotype and its confidence
(posterior-derived GQ as `Phred`; step 13 refines quality later, it does not replace this
field), mint the owned `Genotype` from the winning `GenotypeIdx`, fill `LocusInference` —
expected copies included, because recomputing them downstream from calls gives a different
number. *Depends:* C2. *Source:* spec §2, §9; arch §2.

> **Checkpoint C:** the loop converges, stops, caps, and reports — proven on fixtures at one
> sample and several. Pause for review.

### Milestone D — the `Lg` table and the inert outer skeleton

**D1. The table, built once; the outer rounds, structurally off.**  ☐
`summarise_condition.rs` assembles the whole of spec §2's pseudocode: the `Lg` table built once
per set of slippage numbers by the sibling row builders (contexts per `(read group, candidate)`
looked up from `StratumFits`, hoisted lookup, unhoisted values); the slippage round and
discovery round present as loops whose bodies are unreachable at the defaults (`rounds = 0`,
`Off`), so the default run is one pass through both. The generic path ignores both configs
**structurally** (spec §5.1's closing paragraph), not by half-honouring them. An instrumented
emission-call counter lands here. *Depends:* C3. *Source:* spec §2, §5; arch §1, §5, §6.

**D2. The two cost invariants.**  ☐
Tests: **the emission-call count** equals `candidates × Σ_s (observations in sample s) × builds`
with `builds = 1` at the defaults, independent of pass count — `Σ_s`, not a three-way product,
because a fixture with equal per-sample observation counts is the one shape that hides the bug
(spec §13 test 5, the off half; the re-fit half lands with the bakeoffs plan); **zero
allocations per pass** — the allocation count over a locus is independent of pass count (test
7). **Own commit, do not bundle** — a per-pass rebuild gives identical answers, only slower;
the instrumented counter is the oracle. *Depends:* D1. *Source:* spec §13 tests 5, 7; arch
§Test & bench shape.

> **Checkpoint D:** the table is paid for once and the defaults cost nothing extra. Pause for
> review.

### Milestone E — the input edge: real evidence in, genotypes out

**E1. Evidence shaping.**  ☐
Build `LocusEvidence::Generic` from a `CohortObservation` (the views over `SampleSupport`'s
`(allele, read group)` rows and partials — data-shaping only, no arithmetic) and
`LocusEvidence::Ssr` from the STR generator's `SequenceObservation`s + `SsrDetail`. *Depends:*
A1. *Source:* arch §2, §7; [`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.

**E2. `FrozenParameters` assembly.**  ☐
One constructor per run from the pre-pass's outputs: calibration scales from the accumulator
(scale 1 + `Defaulted` where no rate), contamination views (absent at one sample — absent, not a
fitted zero), per-sample `F` in run order, the seed from `project_spectrum_seed`, the
`StratumFits` borrow, ploidy. **Not class frequencies** — the mixture's second half is per locus
and is built inside the loop, at E2a. *Depends:* A1. *Source:* arch §2; the three parameter
sources' plans.

**E2a. The contaminant frequency, per locus and per sample.**  ☐
Wire the two halves the read likelihood built: `fill_batch_allele_copies` once per locus per
iteration from the loop's current expected copies, then `fill_contaminant_allele_frequencies`
**once per sample**, leaving that sample's own copies out of its own batch, into the
`ContaminationMixture` the row reads. Two buffers on `CallingScratch`, both
`batches × alleles` — the copies, which are per locus, and the frequencies, which are per
sample.

**What this step owes beyond the wiring:** the batching itself. `SequencingBatches`
([`../arch/parameter_prepass_joint_fit.md`](../arch/parameter_prepass_joint_fit.md) §1.6) is
specified and unbuilt, and the loop needs two views of it — `BatchOfEachReadGroup` for the
mixture and `BatchOfEachSample` for the fill. **Whoever builds it owes the rule for a sample
whose libraries ran in different batches**, which the read likelihood deliberately did not
invent (C2's report §6). Until it exists the default — every read group together — is what a run
gets, and that is a complete answer rather than a stub. *Depends:* C1, E2. *Source:*
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §3.6;
[`calling_read_likelihoods.md`](calling_read_likelihoods.md) C2.

**E2b. The run says what contamination it used, per sample.**  ☐
**Spec §3.6 requires it and nothing owned it until now** (added 2026-08-24, on the owner's
instruction, after C1's and C2's reviews both found it homeless): *the run's output must still
carry the fraction used, per sample, because a genotype computed at `c = 0.03` and one at
`c = 0` are otherwise indistinguishable.*

What has to travel, and none of it is a summary: the fraction itself; **whose reads it was
fitted from** (`ContaminationSource` — a library's own reads, or the whole sample's copied onto
it, which are different claims); the two evidence counts beside it, because **a fraction near
zero because nothing could be measured is not the same claim as one measured and found clean**
and only the counts tell them apart; and whether the batching the frequencies were drawn against
was declared or defaulted, which `SequencingBatches::is_default` exists to answer and which the
dense `BatchOfEachReadGroup` the mixture holds cannot (C2's review).

`ContaminationView` already carries the first three. **What is missing is a route from there to
the output**, and a decision this step must take rather than inherit: per sample or per read
group. The fraction is fitted per read group and a sample may hold several, so a per-sample line
either picks one or summarises; §3.6 asks for per sample and the finer grain is what the fit
produces. *Depends:* E2a, E3. *Source:*
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §3.6.

**E3. The integration fixture — ng calls genotypes.**  ☐
End-to-end over a small fixture: reads → merge → E1's shaping → `call_locus` →
`LocusInference` asserted against hand-derived genotypes, on both paths, at one sample and at
several. Provenance and `seed_diversity_unreachable` reach the output. **This is the milestone
where ng calls genotypes.**

**Where the candidates come from, restated 2026-08-25 now that selection exists.** On the
**generic path**, run `select_generic` on the merged locus, so the fixture is reads → merge →
selection → shaping → loop and nothing in it is supplied by hand. **The three joins that step
has to get right are named in
[`../arch/candidate_alleles.md`](../arch/candidate_alleles.md) §5.1** — the prior takes the
*total* allele count and not `alternative_allele_count()`; `LocusSelection::unmatched` is
parallel to the merge's covering samples and not to the run's sample order; and
`genotype_must_be_missing` has no carrier in `SampleGenotypeCall` until this plan adds one.

**The repeat-tract half of this step is blocked, and on the row rather than on the candidates.**
Without the STR read-likelihood row there is no `Lg` for a tract genotype, so no supplied
candidate set rescues it: `censored_emission` is `unimplemented!()`
([`calling_read_likelihoods.md`](calling_read_likelihoods.md) G1) and the row itself is that
plan's H1 and H2. When it lands, this step's tract candidates are still **fixture-supplied**,
because the STR selection path is unwritten — a second and independent gap. Say both in the
test's own doc comment, so a later reader does not read one as the other. *Depends:* D1, E1, E2,
E2a; the STR half additionally on that plan's G1–H2. *Source:* spec §1, §9.

> **Checkpoint E:** genotypes come out of real evidence — over selected candidates on the generic
> path, over supplied ones at a repeat tract. Pause for review.

### Milestone F — the two loop oracles

**F1. The SNP/indel parity oracle.**  ☐
Given the same likelihood table, the same prior parameters and the same candidate set,
production's loop ([`posterior_engine.rs:2733`](../../../../src/var_calling/posterior_engine.rs))
and ng's produce the same genotypes — any difference traced to a decision one of the three
calling documents records, and the trace written into the test. *Depends:* E3. *Source:* spec
§10; §13 test 8.

**F2. The STR differential.**  ☐
On the same likelihood table, run ng's loop under **production's convergence rule and
tolerance** (π at `1e-6`) and require matching genotypes against production's STR loop; then
restore ng's rule (copies over chromosomes at `1e-3`) and **report what moved**. A differential
with a failing state, not parity with an escape clause. *Depends:* E3. *Source:* spec §10.

> **Checkpoint F:** ng's loop is anchored to production on both paths. **Arm A is complete; ng
> calls genotypes.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | compilation against both siblings; ordering-contract assertions |
| B | hand-computed E-step case; **bitwise M-step mutation check on summed copies** |
| C | the flat-pass trap test; one-sample bitwise fixed point; capped-emit flag; **the two-cohort-size division test** |
| D | **instrumented emission-call count** (`candidates × Σ_s obs × 1`); zero-allocation-per-pass |
| E | the integration fixture — hand-derived genotypes from real evidence, both paths; **selected** candidates on the generic path, supplied at a repeat tract |
| F | **production parity (SNP/indel)** and **the STR convergence differential** |

## Out of scope (next plans)

- **Arms B/C, the exhaustive scorer, the local search; the slippage re-fit and discovery
  bodies; every Q1/Q2/Q3 measurement** — [`calling_bakeoffs.md`](calling_bakeoffs.md).
- **Candidate selection (step 6)** — its own plan
  ([`candidate_alleles.md`](candidate_alleles.md)), **built and merged on the generic path**.
  Consume `select_generic`; **do not edit `src/ng/calling/allele_candidates/`**, which is still
  moving through that plan's Milestone D. The repeat-tract half
  ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)) is unwritten, and that is what still
  blocks an end-to-end run at a tract.
- **Wiring `call_locus` into the merge's builder for real runs** — the placement commutes
  (spec §9). With selection merged, this wiring is the only thing between the generic path and
  the GIAB and HG002 regressions the sibling specs name as their definition of done.
- **The deferred loop items** — the shared-prior fast path for large cohorts, the one-sample
  second-pass skip (Q6), threshold/cap tuning from the `passes` distribution (Q4) — spec §11,
  §12.
