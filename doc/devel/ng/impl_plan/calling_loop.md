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

**It is still blocked for repeat tracts, and this paragraph has been wrong about how many gaps
there are twice** — it named two in the original, one after C3b's review of 2026-08-25, and the
count below is the third answer. Each correction was a real finding and each replaced the sentence
rather than appending to it, which is why the history is recorded here rather than in the prose.

**What is not a gap: the row.** The original said the STR read-likelihood row did not exist because
`censored_emission` was `unimplemented!()`. It exists — it is
[`likelihood/ssr.rs`](../../../../src/ng/calling/likelihood/ssr.rs)'s
`genotype_log_likelihood_row` over the shipped `StutterSubstitutionEmission`, landed as
[`calling_read_likelihoods.md`](calling_read_likelihoods.md)'s H1 and H2 and merged; the one
`unimplemented!()` left anywhere under `src/ng/calling/` belongs to a `#[cfg(test)]` oracle that
scores complete observations only.

**What is a gap, as of 2026-08-26, is three things, and E2c closes the middle one.**

1. **The STR candidate path is unwritten**: `allele_candidates/` holds `mod.rs` and `generic.rs`
   only, and [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) is a later session's — so a
   tract's candidates are **fixture-supplied** rather than selected.
2. **Nothing assembled what the row takes.** It takes a scoring context per
   `(read group, candidate)`, and nothing outside `likelihood/ssr.rs`'s own tests had ever built
   one — the gap the owner's ruling of 2026-08-26 turned into **E2c**, which is where the
   sentence *"supplying the candidates is enough"* was retired.
3. **The driver had no route from a tract's evidence to that row** — its emission build and its
   row assembly both took the SNP/indel path's per-sample evidence. **Closed by E3b, 2026-08-27**:
   five places branch, the front-door refusal is gone, and a tract is scored. A *contaminated*
   tract was refused for one step longer, until the third term of its read-likelihood mixture was
   built; **that is E2d, closed the same day**. Nothing on this path is refused now.

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
- **Emission, site filters, phasing** — steps 10–12
  ([`ng_proposal.md`](../spec/ng_proposal.md)).
- **The artifact correction to the site quality, and the output stage that applies it** —
  [`../spec/calling_quality.md`](../spec/calling_quality.md) §3.4, §6. **Its other half is not out:**
  that spec's §3 places the genotype quality and the *uncorrected* site quality inside this plan's
  step C3, because the posterior row and the likelihood table are both gone by the time anything
  downstream runs. C3's entry below says what that adds.
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

**A2. `LocusGenotyper` + `CallingLoopConfig`.**  ✅
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

**B1. The E-step for one sample.**  ✅
Given a filled `Lg` row: build the sample's concentration (`sample_concentration`), its log-prior
row (`GenotypePriorModel`), add, softmax into the posterior row, fold into per-sample expected
copies. Pure, scratch-backed, no allocation. Test: a hand-computed 2-allele diploid case, prior
and likelihood chosen so every intermediate is checkable by hand. *Depends:* A2. *Source:* spec
§2; arch §1.

**B2. The M-step, fixed order.**  ✅
The cohort's expected copies as a sum over samples **in the run's fixed sample order**. Tests
(spec §13 test 2): permuting the samples changes no genotype; **the mutation check is on the
summed expected copies, compared bitwise** — not on the argmax, which will not flip on a
last-bit reorder and proves nothing. **Own commit, do not bundle** — a reordered float sum is
quietly different output at another worker count, never a crash; the bitwise test is the oracle.
*Depends:* B1. *Source:* spec §8; arch §4.

> **Checkpoint B:** one pass is exact, deterministic, and hand-verified. Pause for review.

### Milestone C — the frequency loop

**C1. The flat first pass.**  ✅
The prior-free initialisation (reads only), ported with its reasoning — `EmStepPhase`'s shape
([`posterior_engine.rs:2586`](../../../../src/var_calling/posterior_engine.rs)); it runs at the
start of every outer round, which the skeleton in D preserves. Test (spec §13 test 3): after the
first pass the expected copies reflect the reads and not the seed; **and the trap itself** — at
a locus thin enough that the seed outweighs each sample's reads, a seeded first pass converges
to no-variant where the flat start converges to the variant. *Depends:* B2. *Source:* spec §3.

**C2. Convergence, the cap, and the emitted flag.**  ✅
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

**C3. The final pass.**  ✅
Score every sample once more, take the highest-posterior genotype and its confidence, mint the
owned `Genotype` from the winning `GenotypeIdx`, fill `LocusInference` — expected copies included,
because recomputing them downstream from calls gives a different number.

**The confidence is not this step's arithmetic** *(amended 2026-08-25)*.
[`../spec/calling_quality.md`](../spec/calling_quality.md) §3.1 and §4 own the formula and the
99 cap; this pass calls that function once per sample as it scores it, because the posterior row is
a single reused buffer and computing it afterwards would need the whole table kept. **And two more
outputs land here**, for the same reason — the inputs are gone otherwise: the site quality before
its artifact correction, from the likelihood table once the loop has stopped (§3.2, §5), and the
nine pooled read counts the correction consumes (§3.3). The correction itself and the output stage
that applies it are that plan's, not this one's. *Depends:* C2. *Source:* spec §2, §9; arch §2;
[`../spec/calling_quality.md`](../spec/calling_quality.md) §3.

> **Checkpoint C:** the loop converges, stops, caps, and reports — proven on fixtures at one
> sample and several. Pause for review.

### Milestone D — the `Lg` table and the inert outer skeleton

**D1. The table, built once; the outer rounds, structurally off.**  ✅
`summarise_condition.rs` assembles the whole of spec §2's pseudocode: the `Lg` table built once
per set of slippage numbers by the sibling row builders (contexts per `(read group, candidate)`
looked up from `StratumFits`, hoisted lookup, unhoisted values); the slippage round and
discovery round present as loops whose bodies are unreachable at the defaults (`rounds = 0`,
`Off`), so the default run is one pass through both. The generic path ignores both configs
**structurally** (spec §5.1's closing paragraph), not by half-honouring them. An instrumented
emission-call counter lands here. *Depends:* C3. *Source:* spec §2, §5; arch §1, §5, §6.

**D2. The two cost invariants.**  ✅
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

**E1. Evidence shaping.**  ✅
Build `LocusEvidence::Generic` from a `CohortObservation` (the views over `SampleSupport`'s
`(allele, read group)` rows and partials — data-shaping only, no arithmetic) and
`LocusEvidence::Ssr` from the STR generator's `SequenceObservation`s + `SsrDetail`. *Depends:*
A1. *Source:* arch §2, §7; [`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.

**E2. `FrozenParameters` assembly.**  ✅
One constructor per run from the pre-pass's outputs: calibration scales from the accumulator
(scale 1 + `Defaulted` where no rate), contamination views (absent at one sample — absent, not a
fitted zero), per-sample `F` in run order, the seed from `project_spectrum_seed`, the
`StratumFits` borrow, ploidy. **Not class frequencies** — the mixture's second half is per locus
and is built inside the loop, at E2a. *Depends:* A1. *Source:* arch §2; the three parameter
sources' plans.

**E2a. The contaminant frequency, per locus and per sample.**  ✅

**⚖ Owner's ruling, 2026-08-26, taken when E2 finished and this step's shape turned out to
contradict D2's.** [`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §3.6 says that with
contamination on *"a caller may no longer cache a whole row across iterations"* — because `q(o)`,
the contaminating population's frequency for the allele an observation shows, is the locus's own
number and moves with the loop. D2 pins the opposite: the genotype-likelihood table built **once**
per locus, asserted as `EmissionCost { table_builds: 1 }`.

**The build splits in two.** The **emission** values — one per `(sample, observation, candidate)` —
read no frequency and stay computed once per locus. The per-genotype **row assembly** runs once per
pass, **and only when contamination is on**. So D2's invariant is kept in the form it was really
about — the expensive half, `candidates × Σ_s (observations in sample s) × builds`, is independent
of the pass count — and an uncontaminated run keeps today's behaviour exactly, which is why D2's
existing fixtures pass unchanged.

**This step therefore edits D1's table build and D2's instrument**, which is more than the entry
below describes. What `EmissionCost` counts has to say which half it is counting.
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

**E2c. The repeat tract's scoring parameters — the assembly nobody had written.**  ✅

**⚖ Owner's ruling, 2026-08-26, taken when E2a finished.** The generic path is ready end to end
and the repeat-tract path is not, and what is missing is neither the row nor a parameter but the
**assembly between them**. The row takes a scoring context per `(read group, candidate)`
([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4), and when the ruling was taken
nothing outside `likelihood/ssr.rs`'s own tests had ever built one: `SsrScoringContext::new` had
no production caller, `fill_reachable_lengths` had none, and the outlier weight had no source at
all. **Build it as its own step, before E3** — E3's job is showing that genotypes come out of real
evidence, not inventing the parameters they come out of.

**What it composes, all of it already built.** `StratumFits::at` for the stutter numbers, keyed
by the **candidate's** repeat count and never the tract's (§4.4); `stutter_model_for` for the
conversion into the seven shares; `FrozenParameters::ssr_substitution_rate_at` for the fourth
fitted number, which E2 made reachable; `fill_reachable_lengths` for the length support the
outlier weight is spread over (§4.5). **It has a borrowing shape of its own** — the contexts
borrow the stutter models, so the models and the contexts cannot live in one struct.

**What it owes beyond the wiring: three answers, and only the first is written down.**

- **Where the outlier weight comes from.** §4.5 settles the number — `DEFAULT_OUTLIER_WEIGHT`,
  0.01, *inherited from production and declared inherited*, with no source in the parameters fit.
  What is not settled is whether that inheritance enters the locus's warrant. **It must not**:
  the warrant is per `(read group, candidate)` and this is one run-wide constant, so folding it
  in would make **every** repeat tract's call `Defaulted` and erase the fitted-against-borrowed
  distinction §4.4 says the warrant exists to carry. The same line is already drawn for
  `PART_REPEAT_SHARE_OF_WHOLE`, a placeholder inside every fitted stutter model that no
  provenance mentions.
- **What a candidate whose stratum the fit never reached is scored under.** `NoSlippage`'s own
  documentation says a caller owes an answer — *"a candidate several repeats from its reference
  tract's length can land here on perfectly good data"* — and names four different absences.
  Nothing in the three calling documents rules on it. **Answered** with
  `StutterModel::hipstr_shipped()` and a `Defaulted` warrant, with the two absences that mean
  *the run is not what it claims* counted apart from the two that are ordinary.
- **The same for the substitution rate**, whose emitter records the gap in so many words: *"there
  is no rung below [`FittedHere`] for this parameter … a case the design has not ruled on"*.
  **Answered** with a stated constant, defined as the SNP/indel path's default so that a run
  cannot default its two error parameters to two different guesses.

**What it did not build**, recorded because the boundary was a sentence rather than a discovery:
the **contaminant seed at a tract**. While it was missing, this assembly **refused a run whose fit
found contamination** rather than handing back the two-term form — a mechanism rather than a doc
comment, because the two-term row returns perfectly plausible numbers with the fitted fraction
silently dropped. **E2d built the seed and retired that refusal**; what is left in its place is a
check that the seed and the fractions came from one run. *(When E2c landed the driver also still
turned away every repeat tract, having no route from a tract's evidence to the row; E3b built
it.)*
*Depends:* D1, E1, E2. *Source:*
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4.2, §4.3, §4.4, §4.5;
[`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §4.1, §4.2.

**E2d. The contaminant seed at a repeat tract.**  ✅
**Ran after E3b, not before it**: E2c refused a contaminated run by name, so E3b's tract fixture
did not need this and nothing else reached a tract yet. The third term of §4.5.1's mixture, which
was the one field of the row's locus parameters E2c left empty.

The prior's seed shape is built **per candidate**
(`genotype_prior::seed_ssr::fill_seed_share_per_candidate`) and `c · seed(o)` asks for a
probability per observed **length**; converting the first into the second over
`fill_reachable_lengths`' support is *"the calling loop's job … the only place that holds both
the candidate table and this support"* (`SsrContaminationMixture::contaminant_length_frequencies`).
The conversion itself is settled there — two candidates spelling one length sum into one entry, a
length no candidate reaches gets nothing and its reads fall to the outlier floor, a read that ran
out gets the mass at or above what it witnessed. **What is *not* settled and must not be decided
here is `calling_priors.md` §5's open question 3**, how two candidates spelling one length should
share that length's mass in the *prior*; keying to lengths is what keeps that one question in one
place. *Depends:* E2a, E2c. *Source:*
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4.5.1;
[`../arch/calling_priors.md`](../arch/calling_priors.md) §5.

**⚑ What this step turned up, and it changes a cost claim rather than a number.** A repeat tract's
contaminant term is the **fit's** length spectrum, frozen before calling — §4.5.1 weighed the
cohort's own per-locus frequencies against it and refused them, because contamination must not
move from one pass to the next. So **a contaminated tract's genotype-likelihood table is still
built once**, where a contaminated ordinary site is assembled again at the head of every pass
(§3.6). The driver's per-locus flag is therefore *does this locus's table move as the loop
iterates*, which is contamination **and** the SNP/indel path — not contamination alone. E2a's
entry above describes the ordinary-site half and is unchanged; this is the tract half it did not
have to consider.

**E2b. The run says what contamination it used, per sample.**  ✅
**Runs after E3a, not before it** — its own *Depends* line says so, and it is listed here because
it belongs beside E2a rather than because it comes next.
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

`ContaminationView` already carries the first three. **What was missing is a route from there to
the output**, and a decision this step took rather than inherited: per sample or per read group.
The fraction is fitted per read group and a sample may hold several, so a per-sample line either
picks one or summarises; §3.6 asks for per sample and the finer grain is what the fit produces.

**⚖ Decided: a row is a read group, and it names the sample it belongs to.** Every sample appears,
with each of its read groups under it, which is what §3.6 asks for; nothing is picked and nothing
is averaged, which is what the finer grain is for. A plant sequenced once, in one lane — every
sample of both benchmark cohorts here — gets exactly one row. **A read group is not a library**:
`@RG LB` is a grouping key several read groups can share, so the row names the read group's own
`@RG ID` beside it, and two lanes of one preparation are two rows a reader can tell apart.

**What it carries beyond the four the entry above names**, because two further items were recorded
as owed here and a third arrived with E2d: **the outlier weight, stated once per run as inherited
rather than fitted** (folding it into a tract's per-cell warrant would mark every tract of every
run `Defaulted` and erase the distinction that warrant carries); **how many `(read group,
candidate)` cells a tract defaulted because the fit does not describe this run's read groups**,
counted apart from the ordinary absences because it is the one that means the parameters and the
reads came from different runs; and **whether a tract's contaminant seed was used**, beside the
rung of the tract ladder its length spectrum came from.

**Those three are per locus and the first four are per run, so they travel separately**: a
`RunParameterReport` for what is frozen before calling starts, and `RepeatTractProvenance` on the
locus for what a tract's own parameters rested on. **The latter replaces
`LocusInference::length_spectrum_rung`** — the rung and the counts are one statement about one
tract, and two optional fields with one rule between them could disagree.
[`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §2's sketch of that type was **three
generations behind the shipped one** and has been brought up to it as a transcription, with a note
naming which of C3b, E2e and E2b changed what. *Depends:* E2a, E3a. *Source:*
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §3.6, §4.5.

**E3. The integration fixture — ng calls genotypes.**  **split into E3a ✅ and E3b ✅, 2026-08-26**

**Why it split, and it is a change of character rather than a blockage.** E3a is a **test**: the
generic path needed no library code that did not already exist, and the step is one fixture.
**E3b is implementation** — the driver's emission build, its row assembly and its evidence
accessor all take the SNP/indel path's per-sample evidence unconditionally, so a tract needs
source changes in five places before any fixture can run. Doing both under one step id would have
put a fixture and a subsystem behind one checkbox.

**It carried one open question and no longer does.** Where a tract's prior belief comes from was
undesigned; it is settled in
[`../spec/population_diversity.md`](../spec/population_diversity.md) and built by **E2e**, which now
sits between E2d and E3b.

**E3a. The generic path, end to end.**  ✅
One cohort locus's per-sample observations → the merge's allele unification and read attribution
→ `select_generic` → `shape_generic_locus` → `call_locus` → `LocusInference`, asserted against
hand-derived genotypes at one sample and at several. **This is where ng calls genotypes.**
*Depends:* D1, E1, E2. *Source:* spec §1, §9;
[`../arch/candidate_alleles.md`](../arch/candidate_alleles.md) §5.1.

**Where the candidates come from, restated 2026-08-25 now that selection exists.** On the
**generic path**, run `select_generic` on the merged locus, so the candidates are chosen rather
than supplied. **The three joins that step
has to get right are named in
[`../arch/candidate_alleles.md`](../arch/candidate_alleles.md) §5.1** — the prior takes the
*total* allele count and not `alternative_allele_count()`; `LocusSelection::unmatched` is
parallel to the merge's covering samples and not to the run's sample order; and
`genotype_must_be_missing` needs a carrier in `SampleGenotypeCall` — **added by C3b**, which is
`SampleGenotypeCall::Missing`.

**The tract's candidates are supplied where the generic path's are chosen**, and E3b's fixture
must say which of the two it is so that a later reader does not read one as the other. **What is
unwritten is the STR selection path** ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)) —
not the row, which landed with [`calling_read_likelihoods.md`](calling_read_likelihoods.md)'s H1
and H2, and not its scoring parameters, which landed with E2c.

**What E3a runs and what it supplies, because the difference is the test's whole worth.** It
supplies the per-sample **observations** — what each sample's reads showed, as the locus generator
emits them, since turning aligned reads into observations is step 5's and outside this plan's
Scope — and it supplies the **`ClosedLocus`** the merge is handed, because `LocusCloser`'s
chaining walk is a different subsystem. What it reproduces rather than assumes is that walk's
**keep rule**, so a fixture the real walk would discard as too quiet is refused rather than
scored; that check caught one such fixture on its first run. Everything after the closed locus is
run: the merge's allele unification and read attribution, `select_generic`,
`shape_generic_locus`, and the loop. It lives in `tests/` so that the seams this path names are
`pub`.

**E2e. The repeat tract's prior seed, from the fit rather than from a construction.**  ✅

**⚖ Design settled 2026-08-26 in [`../spec/population_diversity.md`](../spec/population_diversity.md),
which supersedes [`../spec/calling_priors.md`](../spec/calling_priors.md) §5's construction.** A
tract's prior belief about which lengths are plausible is no longer built as a geometric decay away
from the cohort's commonest length and scaled to reproduce a measured diversity. It is read from
what the joint repeat fit already produces per stratum: a **length spectrum** (the shape) and a
**concentration** (the strength), both fitted from that stratum's own tracts and already run end to
end on both benchmark cohorts.

**Why this is a step of this plan and not only of the prior's.** It is the one thing between the
loop and a called tract, and it touches three modules — the fit's seam, the prior's seed builder,
and the run's frozen parameters. **It edits `genotype_prior/seed_ssr.rs`, which is
[`calling_prior.md`](calling_prior.md)'s module**; recorded as a deviation rather than done quietly,
because that plan's E1 shipped the function this replaces.

**Three pieces:**

1. **Carry the two values across the seam.** `StratumFits`
   ([`joint/stratum_fits.rs`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)) gathers
   each stratum's slippage numbers and their provenance and **drops the length spectrum and the
   concentration** that `fit_strata` produced. They are carried, keyed the same way, with the rung
   they came from beside them.
2. **Rebuild the seed on them.** `fill_ssr_seed` takes the fitted pair and the locus's candidate
   lengths, and maps the spectrum — indexed by offset from the **reference** tract length — onto
   them. **Two inputs disappear**: the cohort's modal repeat count at the tract, which had no source
   because repeat-tract selection is unwritten, and the run-wide repeat gene diversity, which
   nothing emits. `SsrSeedOutcome::DiversityUnreachable` goes with them — the failure that fires at
   *every* tract at one outbred sample exists only because a constructed shape had to be scaled to a
   measurement.
3. **The three-rung fallback**, spec §4.4: the stratum's own fit; failing that its motif period's
   pooled tracts; failing that a flat spectrum at a stated concentration. Each marked, each
   reported. **Measured, and it is why the rung is a pooled fit rather than a curve**: the strata
   with no fit of their own hold about 2% of HG002's tracts and at most 7% of tomato's.

**What it does not do:** the ordinary-site side, which is a separate and smaller job and blocks
nothing here. **Both its numbers are already fitted** — the joint fit emits the population's
allele-frequency density and the heterozygosity read off it
([`joint/fit.rs:207`](../../../../src/ng/parameter_estimation/joint/fit.rs), `:223`) — and what
stands between them and the caller is a **representation change**: wrap the heterozygosity in its
newtype, and project the density into the `2N + 1` allele-count classes the seed projection takes.
**Home:** a step of its own, spec §3.2. *Depends:* E2c. *Source:*
[`../spec/population_diversity.md`](../spec/population_diversity.md) §4, §5.

**E2f. The ordinary-site prior's seed, from the fit.**  ✅

**Both numbers are already fitted and neither reaches the caller.** `JointFit`
([`joint/fit.rs:198`](../../../../src/ng/parameter_estimation/joint/fit.rs)) carries the
population's allele-frequency density (`:207`) — a Beta over the segregating positions plus a point
mass at each end — and the expected heterozygosity read off it (`:223`). The caller's seed
projection takes an `ExpectedHeterozygosity` and a `FittedSpectrum` of `2N + 1` allele-count class
weights, and its own seam is already built: `RunParameters::project_seed`
([`run_parameters.rs:97`](../../../../src/ng/calling/run_parameters.rs)) takes both as arguments and
nothing supplies them. **So this is a representation change, not an estimator.**

**Two pieces:**

1. **Wrap the heterozygosity** in its newtype, whose constructor rejects a value outside `[0, 1]`.
2. **Project the density into allele-count classes** at the panel's own size — the Beta evaluated
   into `2N + 1` classes with the two point masses at the ends — plus the two bookkeeping numbers
   `FittedSpectrum` carries.

**⚑ ANSWERED elsewhere, 2026-08-26 — this step has no panel-size floor to decide.** The entry used
to say that without a floor *"the ordinary-site ladder's top two rungs are separated by nothing"*
and point at spec §9's question 3. Both halves are retired: the top two rungs are now the two ends
of one ramp, a run slides between them at a weight that rises with the panel, and §9's question 3
is answered *no floor* on the branch `ng-seed-shrinkage`. See
[`../spec/ordinary_site_seed.md`](../spec/ordinary_site_seed.md) §4 and
[`../../reports/implementations/ng_seed_shrinkage_2026-08-26.md`](../../reports/implementations/ng_seed_shrinkage_2026-08-26.md).

**What this buys, and it is worth stating because the step is small:** the SNP/indel prior moves off
`ExpectedHeterozygosity::SPECIES_FALLBACK` — a species-range constant — onto this cohort's own
measurement, and gains a fitted shape wherever the panel supports one. *Depends:* E2. *Source:*
[`../spec/population_diversity.md`](../spec/population_diversity.md) §3.

**E3b. The repeat-tract path, end to end.**  ✅
The same fixture at a tract, with the candidates and their repeat counts **fixture-supplied**
rather than selected — the STR selection path is unwritten
([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)), so the test's own doc comment must say
which of the two it is, and a later reader must not read a supplied candidate set as a selected
one.

**The implementation, which needs no ruling.** Five places in
`inference/summarise_condition.rs` take the generic path's per-sample evidence unconditionally and
each needs a tract branch: `generic_evidence_of`, `build_locus_emissions`,
`assemble_genotype_likelihood_table`, `weakest_warrant_at_the_locus` and the locus concentration
fill; `call_locus` refuses a tract in front of all of them. The tract's own row and its scoring
parameters are built (E2c), and with no contamination the row reads no frequency, so the table is
assembled once per locus and D2's invariant holds unchanged.

**The one per-locus input a tract needs beyond its candidates**, fixture-supplied here because the
STR selection path is unwritten: each candidate's **repeat count**, which is not derivable from its
bases — an interrupted tract holds fewer whole repeats than its length suggests. *(The cohort's
modal repeat count at the tract was a second such input until E2e retired it: the fitted spectrum is
indexed by offset from the reference tract length, which every locus already knows.)*

**⚖ The open question this entry carried is closed**, by
[`../spec/population_diversity.md`](../spec/population_diversity.md) and step E2e above: the tract's
prior is seeded from the per-stratum length spectrum and concentration the fit already produces,
carried on the run's frozen parameters, present or absent as a whole. Neither the cohort-wide repeat
diversity nor the decay per repeat is needed. **What E3b still owes is the driver's tract branch and
the fixture.**

*Depends:* E3a, E2e. *Source:* spec §1, §5, §9;
[`../arch/calling_priors.md`](../arch/calling_priors.md) §5.

**⚖ Owner's ruling, 2026-08-27, on the one point where this step and a spec section disagreed.**
[`../spec/population_diversity.md`](../spec/population_diversity.md) §5 wanted a tract in a run
carrying no repeat-tract parameters **refused by name**; §4.4 wants the tract ladder to always
answer. The two meet only at a run whose fit produced no length spectrum anywhere. **The tract is
called**, and its record carries the bottom rung — *"refusing turns a whole class of runs into a
hard failure for a condition the output already states"*. §5 and §6 of that spec now record the
ruling; nothing in the code changed, because this is what E3b built.

> **Checkpoint E:** genotypes come out of real evidence — over selected candidates on the generic
> path (E3a, done), over supplied ones at a repeat tract (E3b, done). Pause for review.

### Milestone F — the two loop oracles

**F1. The SNP/indel parity oracle.**  ☐
Given the same likelihood table, the same prior parameters and the same candidate set,
production's loop ([`posterior_engine.rs:2733`](../../../../src/var_calling/posterior_engine.rs))
and ng's produce the same genotypes — any difference traced to a decision one of the three
calling documents records, and the trace written into the test. *Depends:* E3a. *Source:* spec
§10; §13 test 8.

**F2. The STR differential.**  ☐
On the same likelihood table, run ng's loop under **production's convergence rule and
tolerance** (π at `1e-6`) and require matching genotypes against production's STR loop; then
restore ng's rule (copies over chromosomes at `1e-3`) and **report what moved**. A differential
with a failing state, not parity with an escape clause. *Depends:* E3b. *Source:* spec §10.

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
| E | the integration fixture — hand-derived genotypes from real evidence; **selected** candidates on the generic path (E3a), supplied at a repeat tract (E3b); and, behind the tract's half of it, an assembly whose every lookup changes a row when it is dropped |
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
