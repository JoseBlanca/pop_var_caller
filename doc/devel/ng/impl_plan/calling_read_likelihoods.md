# ng read likelihoods (step 7) — implementation plan

**Status:** draft, 2026-08-21. The build order for **step 7 whole**: the
`calling/likelihood/` module — the evidence views, the SNP/indel closed form with its
contamination mixture and partial-read rule, the stutter distribution's three recorded changes,
the STR emission seam with Model A and the censored term, and the STR row with its three-term
mixture. Design is settled in [`read_likelihoods.md`](../spec/read_likelihoods.md) (spec) and
[`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §1–§4 (types & interfaces). This
plan turns that design into build order; it is **not** a place for new design — the open
questions (Q1, Q7) are the spec's, and neither blocks a step below.

*(One plan, not a per-path split: the shared milestones (A) are small and the two paths land as
disjoint milestone runs (B–D generic, E–H STR), so a reader can execute one path's run without
the other. If execution proves the file unwieldy, split it exactly as
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) /
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) are split — the milestone boundaries below
are the cut line.)*

**Where this sits.** Six plans build calling:
[`calling_prerequisites`](calling_prerequisites.md) ∥
[`calling_foundations`](calling_foundations.md) → [`calling_prior`](calling_prior.md) ∥
`calling_read_likelihoods` → [`calling_loop`](calling_loop.md) →
[`calling_bakeoffs`](calling_bakeoffs.md). **This plan needs prerequisites items 1, 2, 3 and 5** —
the merge's read-group axis and kept partials are what the generic evidence views read; the
calibration accumulator and the `StratumFits` gather are what its parameter views consume. The
prior plan needs none of that, which is why it starts earlier; the two still run in parallel once
this one starts.

**Item 4 came here instead of arriving from there (owner, 2026-08-24).** The contamination
mixture's second half was to be three allele-class frequencies emitted by the parameter pre-pass;
it is now the frequency of the allele an observation shows, at the locus being called, over the
samples in this sample's sequencing batch — **recomputed every iteration**, because it is a
property of the locus and the loop already estimates it. **So this plan builds it**, in Milestone C,
and the pre-pass owes nothing (spec §3.6, corrected the same day;
[`calling_prerequisites.md`](calling_prerequisites.md) Milestone E records why it could never have
been the pre-pass's).

---

## Scope

**In:** `src/ng/calling/likelihood/` — `mod.rs` (evidence views, the row contract, the parameter
tiers, shared floors, scratch shapes), `generic.rs` (the SNP/indel closed form + contamination
mixture + the partial rule), `ssr.rs` (the STR row), `ssr_emission.rs` (`SsrEmissionModel`,
`StutterSubstitutionEmission`, `ClassicEmissionOracle`); the three recorded changes to
[`alignment/stutter.rs`](../../../../src/ng/alignment/stutter.rs) (ng code, not frozen
production); the two doc repointings spec §7 asks for.

**Out (later plans or upstream):**

- **When rows are built, cached and reused; the `Lg` table; context assembly per locus** — the
  loop's ([`calling_loop.md`](calling_loop.md)).
- **The change measurements** — the dropped coefficient, the STR contamination on/off sweep, the
  STR-vs-generic zero-slippage comparison, the per-locus re-fit study (spec §12 items 14–19).
  All need genotypes end-to-end; recorded in [`calling_bakeoffs.md`](calling_bakeoffs.md)'s
  out-of-scope as unscheduled measurement runs.
- **`q_sum_other`'s producer.** The pooled leftover the formula keeps is created by candidate
  *selection*, which has no spec — "whoever specifies selection owes the pool"
  ([`calling_em_loop.md`](../spec/calling_em_loop.md) §4). The row takes it as an input
  (`unmatched_q_sum`); tests supply it from fixtures; **no producer is invented here.**
- **A part-repeat estimator, untying the one-step shares, the purity adjustment, per-base STR
  qualities, an overdispersion device** — deferred with homes (spec §10).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** Each path's scoring mathematics is built and
  proven against an oracle on hand-built evidence before anything touches the merge's real
  output; the generic path's run (B–D) is ordered plain form → mixture → partials so each step's
  oracle is the previous step.
- **Reuse over rewrite.** The closed form's shape, the mixture, the placement enumeration, the
  model seam and the comparator are **ports** (spec §9's map); the stutter distribution is
  **reused from `alignment/stutter.rs`**, renamed — one built implementation with two consumers,
  never a second spelling. The substitution comparison is **composed** from `FlatEmission`
  ([`alignment/emission.rs:250`](../../../../src/ng/alignment/emission.rs)), never re-implemented.
  HipSTR's GPL tree is **not a source to implement from**; the operative source is spec §4.2's
  own full statement of the distribution (the licence rule, spec §4.2).
- **Verify against ground truth.** The generic row against production's
  `standard_log_likelihood`
  ([`per_group_merger.rs:1948`](../../../../src/var_calling/per_group_merger.rs)) with the two
  recorded changes reconciled term-by-term; the STR emission against `ClassicEmissionOracle` —
  an independent implementation, ported test-only exactly as production keeps it.
- **Isolate the silent steps.** The `m(a, g)` divisor, an inverted one-step share, the
  aggregation identity and the outlier spread's cohort independence all fail as quietly-wrong
  genotypes; each lands as its own commit with its oracle named, marked below.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place)

- **Foundations merged:** `GenotypeTableView`, `CandidateAlleles`, `AlleleId`
  ([`calling_foundations.md`](calling_foundations.md)).
- **Prerequisites merged, items 1, 2, 3 and 5**
  ([`calling_prerequisites.md`](calling_prerequisites.md)): the `(allele, read group)` support rows
  and kept partials (Milestones B–C there), the calibration accumulator (D — on the histogram
  route; the census route's waits on the comparison between the two routes), the `StratumFits`
  gather (F). **Item 4 is not owed by that plan** — see the note above.
- The STR evidence needs no merge work: `SequenceObservation`
  ([`locus_generation/mod.rs:295`](../../../../src/ng/locus_generation/mod.rs)) already keys on
  `(bases, witness, read_group)`; `complete_observations()`
  ([`:134`](../../../../src/ng/locus_generation/mod.rs)) is where the split is spelled there.
  *(Corrected 2026-08-24 at A1: this said that method "stays the only unguarded access", and it is
  neither a guard — the field it reads is `pub` — nor now the only one, since `SsrSampleEvidence`
  holds a bare slice and spells the split again. Arch §2.2 carries the full correction.)*
- The pieces to reuse or port: `StutterModel`/`StutterRates`/`probability`
  ([`alignment/stutter.rs:147`](../../../../src/ng/alignment/stutter.rs),
  [`:82`](../../../../src/ng/alignment/stutter.rs),
  [`:300`](../../../../src/ng/alignment/stutter.rs)), `MAX_SLIP = 10`
  ([`:63`](../../../../src/ng/alignment/stutter.rs)), the geometric clamps
  ([`:74`](../../../../src/ng/alignment/stutter.rs), made `pub` at A2 so the likelihood names them
  rather than spelling a second copy); production's seam and comparator
  ([`read_model/mod.rs:63`](../../../../src/ssr/cohort/read_model/mod.rs),
  [`classic.rs`](../../../../src/ssr/cohort/read_model/classic.rs)); the placement enumeration
  ([`ssr/cohort/stutter.rs`](../../../../src/ssr/cohort/stutter.rs)); `MIN_BASE_ERROR`
  ([`contamination_estimation.rs:1449`](../../../../src/var_calling/contamination_estimation.rs)).

## Worktree, branch, merge

- **Worktree** `../pop_var_caller-calling-read-likelihoods`, **branch**
  `ng-calling-read-likelihoods`, from `main` **after both phase-1 branches
  (`ng-calling-foundations`, `ng-calling-prerequisites`) have merged**.
- **Runs in parallel with** `ng-calling-prior`. Conflict surface: `src/ng/calling/mod.rs` — one
  `pub mod likelihood;` line plus re-exports, placed alphabetically. This branch adds nothing to
  `types.rs`; its floors are named constants in `likelihood/mod.rs` and the stutter edits stay in
  `alignment/stutter.rs`, which no other calling branch touches.
- **Merge order back: the prior merges first, this branch second**, resolving any adjacent-line
  `mod.rs` conflict. The loop plan branches only after both are in.

---

## The steps

### Milestone A — module scaffold + the shared vocabulary (types, no logic)

**A1. Scaffold + evidence views.**  ✅ *(shipped with the merge-to-candidate mapping as an argument
rather than an assumption — `GenericObservation::fill_from_supported_alleles` takes selection's
`&[Option<AlleleId>]` and returns the quality of the rows it dropped, because the merge's allele
index and `AlleleId` are two numberings and a prune renumbers between them; arch §2.1 carries the
correction)*
`calling/likelihood/mod.rs` wired into `calling/mod.rs`. `GenericSampleEvidence`
(`supported: &[GenericObservation]`, `unmatched_q_sum`, `partials`) and `GenericObservation`
(`allele`, `read_group`, `num_reads`, `q_sum`) as **views over the merge's rows** — one entry per
`(allele, read group)`, which prerequisites Milestone B made real; `PartialObservation` over
Milestone C's kept rows, bases + witnessed positions intact. `SsrSampleEvidence` — a slice of
`SequenceObservation` plus the locus's `SsrDetail`
([`locus_generation/mod.rs:438`](../../../../src/ng/locus_generation/mod.rs)). *Source:* arch
§2.1, §2.2; spec §1.4, §2.3.

**A2. Parameter views, floors, contract.**  ✅ *(shipped with three departures, all recorded in the
step's report: the row's scratch is **two** types, because the evidence view borrows the staging
buffer and a row taking `&mut` the same object cannot be called — so `GenericEvidenceBuffer` holds
the evidence and the generic row's own scratch arrives at D1 with the buffers it needs; the
emission cache takes the evidence rather than an observation count, because "how many observations"
had two reachable readings and one of them is silent; and the depth-cap question spec §3.2 hands
to this step is decided — **the denominator stays unthinned** — with the owner's ruling invited)*
`ReadGroupCalibration { scale, provenance }` (scale 1.0 + `Defaulted` where no rate was emitted —
visible, never silent) and `ContaminationView` (fraction, `markers_with_reads`, `reads_on_markers` —
"measured clean" vs "unmeasurable" told apart by the counts; *the plan said `reads_at_markers`,
which is not the field's name — corrected 2026-08-24*). **`ContaminationView` carries no
allele-class frequencies**: the mixture's second half is per locus and per iteration, so it is read
where the allele frequency is rather than frozen on this view (C2 below).

**The scale is `fitted rate ÷ exp(Σ q_sum / Σ num_obs)`, a geometric mean** — the pre-pass's
accumulator supplies the two sums
([`generic/calibration.rs`](../../../../src/ng/parameter_estimation/generic/calibration.rs)), and
spec §3.2 records why that average and not the arithmetic one. Import `MIN_BASE_ERROR = 1e-12` and
the geometric clamps as named constants with their reasons. The row contract in `mod.rs`'s docs: pure function, fills caller scratch, empty evidence
row = all zeros (the prior decides, no branch), mis-shaped input = assertion held in release
([`per_group_merger.rs:1963`](../../../../src/var_calling/per_group_merger.rs) is the precedent).
Scratch shells `GenericRowScratch`, `SsrRowScratch<S>`. *Depends:* A1. *Source:* arch §1.1, §2.3;
spec §3.2, §3.6, §8.

> **Checkpoint A — reached 2026-08-24.** The vocabulary compiles; the tier table is documented on
> the types. *Its three tiers are **frozen / per-call / per-iteration**, not "frozen / per-call /
> invisible" as this line said: spec §3.6's correction of 2026-08-24 gives the third tier a reader
> in this module, and the spec's own §6.1 is corrected to match at A2.* Pause for review.

### Milestone B — the SNP/indel closed form

**B1. `error_spread_divisors` — `m(a, g)`.**  ✅ *(shipped as `fill_error_spread_divisors` plus
`DivisorTable`, which carries the stride: a bare `(values, allele_count)` pair cannot check that
the count is the stride the buffer really has, and reading a three-allele table at a stride of two
returns a real divisor from the wrong row on six of twelve lookups with nothing to panic about —
which is this step's own named failure shape arriving through the accessor written to prevent it.
Verified against an independent oracle over 1,758,811 cells at ploidy 1 to 4: no disagreements.)*
`3.0` where the observation differs from every allele the genotype carries by a substitution at
exactly one position, `1.0` otherwise — a property of the allele pair, computed once per
`(allele, genotype)` over the projected sequences the merge unified. **Own commit, do not
bundle** — a wrong divisor is `log 3` (4.8 Phred) per wrong read in the wrong direction, and
nothing crashes; the oracle is a hand-built fixture per class (one-substitution, multi-position,
indel). *Depends:* A2. *Source:* spec §3.5; arch §3.

**B2. `genotype_log_likelihood_row`, plain form.**  ✅ *(shipped with the error-spread table
storing `log m` — the decision this step owed — and with the production differential reconciling
**three** differences rather than two: production has no calibration, so ng's `n·log scale` comes
back out too, and at a defaulted calibration that term is exactly zero, which is how a mutation
deleting it outright survived every test at a cost of 1,620 Phred. Two specification claims were
found false and corrected: the aggregation identity is not bitwise and neither is order
independence — §12 tests 8 and 9, §2.3, and the architecture's own row contract.)*
Spec §3.3 exactly, with an empty contamination slice: explained reads charged `n·log(k_a/P)`;
unexplained charged `q_sum + n·(log scale − log m)`; `unmatched_q_sum` added as the
genotype-independent constant (kept for emission). Tests: **the aggregation identity, bit for
bit** — the likelihood from a list of individual reads and from their merged aggregate agree with
no round trip through probability space (spec §12 test 9 — the test that would have caught the
geometric-mean substitution, and the reason the formula has this shape); order independence
(test 8); empty row = zeros; a hand-computed biallelic diploid case. **The production
differential:** on the same inputs, ng's row plus the multinomial coefficient (computed in the
test in closed form) minus the `÷3` effect equals `standard_log_likelihood`
([`per_group_merger.rs:1948`](../../../../src/var_calling/per_group_merger.rs)) to
floating-point tolerance — every difference attributed to the two recorded changes (spec §3.4,
§3.5), none unexplained. **Own commit, do not bundle.** *Depends:* B1. *Source:* spec §3.3; arch
§3.

**Three things B1's review put on this step** (2026-08-24), the first of them a gate:

- **The differential must take the `÷3` effect from `DivisorTable`, not from a literal
  `n_alt · ln 3`.** Written with the literal it passes with B1 deleted and the divisor hardcoded —
  the same shape as a test B1 shipped and had to repair, which computed `3.0_f64.ln()` and never
  touched the table. More generally: **at least one B2 test must obtain its divisors by calling
  `fill_error_spread_divisors` on a real candidate table**, so that deleting B1 is a compile error
  rather than a quiet subtraction. Every step of this plan so far reverts green, and this is where
  that stops.
- **Decide whether the table stores `m` or `log m`.** The row charges `log m` once per
  `(observation, genotype)`, so as it stands B2 calls `.ln()` in the inner loop — measured at 1.392
  ns a term against 0.553 ns if the table held the logarithm, about 26 s single-threaded over a
  high-depth sample. It has to be decided here because storing the logarithm makes *divisor* the
  wrong word (arch §3 carries the `OPEN:`).
- **Widen `standard_log_likelihood` to `pub(crate)`** ([`per_group_merger.rs:1948`]
  (../../../../src/var_calling/per_group_merger.rs)) — it is a private `fn` and the differential
  cannot reach it. Visibility only, which is the one production change the freeze allows.

> **Checkpoint B — reached 2026-08-24.** The closed form matches production term-for-term once the
> **three** recorded changes are reconciled — the dropped coefficient, the error spread, and ng's
> calibration scale, which production has no counterpart for and which this line used to leave out.
> *Aggregation is exact in the model and not in the arithmetic: a relative 2 × 10⁻¹⁴, spec §2.3.*
> Pause for review.

### Milestone C — the generic contamination mixture

**C1. The mixture, no `c = 0` branch.**  ✅ *(shipped as the row's one path, with
`ContaminationMixture` holding both halves rather than the sketch's bare
`&[ContaminationView]` — the fraction and the frequency sit in different tiers, and one
construction is where they can be checked against each other. The `c = 0` agreement is a
relative 7.3 × 10⁻¹⁵ over 4,440 comparisons, and the sweep kills both defects this step names:
production's extra `(1 − ε)` factor disagrees on 3,172 of them, its allele-count divisor on
2,336. **One A2 decision had to be reversed to get here** — `calibrated_error`'s ceiling is a
non-linear function of a per-read quality, which spec §2.3 forbids outright, so what the row
charges floors and does not cap; on the aggregation fixture the cap would have moved the answer
69 nats where the property is pinned to a relative 2 × 10⁻¹⁴. Two open questions for the owner
are recorded in the step's report: whether the capped reading survives at all, and whether the
spread table should now store `m`.)*
Spec §3.6: `n_o · log[(1 − c)·own(o|g) + c·q(o)]`, evaluated in probability space with one
logarithm, `q(o)` taken as a parameter for now so this step is about the mixture and nothing else.
**There is no `c == 0` branch** — the two forms agree to a few ulp and the tolerance is a named
test constant (spec §12 test 11); the test also fails the moment anyone reintroduces production's
extra `(1 − ε)` factor or its allele-count divisor into `own`. A second test hand-computes one
contaminated case (`c = 0.03`, contaminant frequency 1 in 1,000). *Depends:* B2. *Source:* spec
§3.6; arch §3.

**C2. `q(o)` — the contaminating population's frequency for the allele the observation shows.**  ✅
*(shipped as `fill_contaminant_allele_frequencies` — a batch's copies summed over its samples and
divided by their total — plus a batch axis on the mixture's frequency table, read through the
observation's own read group. **The row needs no sample identity for this**, which is what made it
fit: batches are over read groups, and every observation already carries one. `SequencingBatches`
is still unbuilt and is not built here; the mixture takes the batching as a
`BatchOfEachReadGroup`, which that type produces trivially when it lands. **Two things are owed and named in the step's
report:** which batch a sample belongs to when its libraries ran in different ones — that rule
belongs with `SequencingBatches`, so the producer takes the sample's batch as an argument — and
whether the never-seen floor should stay defensive at `1e-12` or become a statistical pseudocount
— **settled 2026-08-24: keep it very low**, since the statistical reading would say *this read
might well be contamination* at every candidate the cohort is thin on. **The review found a third
and the owner settled it the same day:** the frequency was summed over the batch's samples
*including the one being scored*, so a sample alone in a batch explained its own alternative reads
as its own contaminant. It now leaves itself out, as the genotype prior's concentration already
did — `fill_batch_allele_copies` sums per locus, `fill_contaminant_allele_frequencies` subtracts
and normalises per sample. **A sample alone in its batch therefore gets the reference**, which is
the conservative answer: a library with no neighbours has no contaminating population.
[`calling_loop.md`](calling_loop.md) E2a and E2b own what is left — the batching itself, and the
run reporting the fraction it used.)*
**New here on 2026-08-24 (owner); this used to arrive from the pre-pass as three allele-class
frequencies and does not.** `q(o)` is the frequency of the observation's own allele at the locus
being called, over the samples in this sample's **sequencing batch** — the grain
[`parameter_prepass_joint_fit.md`](../arch/parameter_prepass_joint_fit.md) §1.6 already carries,
whose default is one batch holding the whole run, so a run declaring no batches gets the cohort
frequency. **Recomputed every iteration** from the loop's current estimate, which is the same
number the genotype prior reads, so this is a lookup and not a fit.

Tests: the frequency of a candidate the batch never shows is the floor rather than zero (a
contaminant read of an allele nobody carries must not make the term collapse); two batches with
different frequencies at one locus give one sample's reads different `q` from another's; and at one
batch the answer equals the cohort frequency, which is what makes the default lose nothing.
**Own commit, do not bundle** — a wrong `q` is a genotype quietly pulled toward the contaminant's
allele and nothing crashes. *Depends:* C1. *Source:* spec §3.6;
[`calling_prerequisites.md`](calling_prerequisites.md) Milestone E for why it is here.

**The one consequence to carry into the loop plan:** the contamination *fraction* is frozen before
the loop and `q(o)` is not, so the two halves of the mixture sit in different tiers. Spec §6.1's
first tier holds the fraction only.

### Milestone D — partial observations on the generic path

**D1. The compatibility rule.**  ✅ *(shipped, and **the rule turned out to need no positional
restriction at all** — the correction that made the step small. An allele is the whole locus as a
carrier has it, not the reference with gaps, so a read from a carrier shows the start or the end
of that carrier's own sequence: flush left means its bases are a **prefix** of the allele, flush
right a **suffix**, and `WitnessedLocusPositions`' two predicates are documented as exactly those
constraints. **So there is no gather and no buffer sized by the widest witness** — this step's
scratch is the compatibility cache alone. A witness with a hole reaching both borders splits into
a prefix and a suffix at an unknown point, which is checked without trying every split; a witness
flush at neither border is anchored to nothing and constrains nothing, rather than being matched
by content somewhere inside the allele, which would move a genotype on a coincidence. The row
gained `&CandidateAlleles` and the scratch, and its two per-read-group parameters were paired into
`ReadGroupParameters` — which is where their read-group counts are now checked against each
other. **The first version shipped a defect its own review caught and this note keeps:** the rule
branched on the two borders alone, so a witness with a hole that reached only one border was
tested as though its bases were contiguous, and the verdicts **inverted** — the allele agreeing
with the read at every position it saw was charged 14 nats and the one disagreeing was charged
nothing. The rule now decides on the run count too, and *says nothing* wherever a run is anchored
to neither border, which is what those shapes genuinely imply rather than a fallback.)*
An allele is compatible with a partial when its projection **restricted to the positions the read
witnessed** equals the partial's bases; a compatible partial contributes `Σ k_a/P` over the
genotype's compatible alleles; compatible with none → charged as an error with `m = 1`. Exactly
aggregable by construction (the witnessed positions are part of the observation's identity).

**Those positions are a *set of runs with holes*, not one run** (spec §5.3, corrected 2026-08-24),
and this step owns the two consequences. The restricted projection is a **gather**, so it needs a
buffer sized by the widest witness — the generic row's own scratch, which A2 deliberately did not
invent and which lands here. And the witness counts *locus positions* while the partial's bases are
what the read showed over them, so **their lengths are not interchangeable** and neither may index
the other. This step also needs a compatibility cache per `(partial, allele)`, for the reason the
STR path caches emissions: the verdict is read by every genotype.

Tests: a partial compatible with both of a diploid heterozygote's alleles contributes 1 — no
information, correctly; the no-compatible error charge; verdicts identical for pooled reads; **and a
witness with a hole in it, against an allele whose bases differ only inside the hole** — the case a
contiguous-range implementation gets wrong and no single-run fixture can see. *Depends:* B2;
prerequisites C. *Source:* spec §5.3; arch §3.

> **Checkpoint C/D:** the generic path is complete — plain, contaminated, and censored evidence
> all scored, each against a hand-computed or production oracle. Pause for review.

### Milestone E — the stutter distribution's three changes (`alignment/stutter.rs`)

The distribution is **reused, not duplicated**; the file already records the follow-ups its port
deferred. All three changes are in ng code.

**E1. Rename to the spec's vocabulary.**  ☐
`StutterRates`/`StutterModel` fields renamed — `whole_repeat_longer_share`,
`part_repeat_one_step_share`, … — with HipSTR's names kept in doc comments for whoever ports
alongside; *in frame / out of frame* is banned vocabulary (spec §1.3) and the fields currently
carry it (`in_up`, `out_geom`). Mechanical; existing tests green unchanged. Alongside, the two
doc repointings spec §7 records: [`alignment.md`](../spec/alignment.md) §5.2 repointed at the
spec's §4.2 as the distribution's owner, its *in frame / out of frame* wording moved to §1.3's.
*Source:* spec §1.3, §4.2, §7; arch §4.2.

**E2. Two named cutoffs + the reported truncation.**  ☐
`MAX_WHOLE_REPEAT_SLIP = 10` (repeats) and `MAX_PART_REPEAT_SLIP = 10` (base pairs) replace the
single `MAX_SLIP` ([`stutter.rs:63`](../../../../src/ng/alignment/stutter.rs)) — both inherited
from production's provisional 10 and declared inherited. The mass the cutoffs discard is
**computed and reported per candidate** (feeds `SsrScoringContext.truncated_mass_lost`). Test:
the reported loss equals one minus the truncated sum, to floating-point tolerance, across
candidate lengths, periods, and one-step shares over the clamped range — **the test pins that
the loss is computed and surfaced, not that it is small** (spec §12 test 5; its size runs from
2 in a million to 2 in a thousand). **Own commit, do not bundle** — an unreported loss compares
candidates on different scales silently. *Depends:* E1. *Source:* spec §4.2; arch §4.2.

**E3. `stutter_rates_for(&Slippage)` + the sums-to-one tripwire.**  ☐
Seven shares from the fit's three numbers, the placeholders named as placeholders:
`PART_REPEAT_SHARE_OF_WHOLE = 0.05` (production's `OUT_FRAME_REL`) and the two one-step shares
tied to one value — both awaiting an owner (spec §10). Test: **the distribution sums to one**
over the full untruncated support, periods 1–6, direction splits symmetric to 5:1 (spec §12
test 4 — the test that catches a one-step share read as its complement, a mis-set same-length
share, and an off-by-one re-indexing, all three silent). **Own commit, do not bundle.**
*Depends:* E1. *Source:* spec §4.2; arch §4.2.

> **Checkpoint E:** the distribution carries the spec's names, two cutoffs, a reported loss, and
> a proven total. Pause for review.

### Milestone F — the STR emission seam

**F1. The seam types.**  ☐
`ssr_emission.rs`: `SsrScoringContext` (the tier-two seam — every number arrives per call, none
read from global state; carries `stutter`, `substitution_rate` — the fitted per-stratum rate,
**never the SNP ε and never `q_sum`**, spec's closed Q6 — `truncated_mass_lost`,
`weakest_provenance`), `SsrCandidate` (bases + the repeat count that keys the stratum lookup —
the **candidate's** stratum, not the reference's, so contexts are per
`(read group, candidate)` and never hoisted out of the candidate loop), and the
`SsrEmissionModel` trait with `emission` and `censored_emission`. *Depends:* A2. *Source:* arch
§4.1; spec §4.3, §4.4.

**F2. `StutterSubstitutionEmission` — Model A.**  ☐
The stutter factor (E's distribution, via the reused `StutterModel`) times the substitution
factor — **composed** from `FlatEmission`
([`alignment/emission.rs:250`](../../../../src/ng/alignment/emission.rs)) under the fitted rate,
over sequences the stutter factor has already made equal-length. Placement enumeration ported
with production's split, stated: whole-repeat slips enumerate placements with equal weight,
part-repeat changes resize at the tract's end in a single placement. Unreachable slips score
zero. Tests: the three ported property tests — same-length dominance (spec §12 test 1),
direction and size ordered as fitted (test 2), a whole repeat beats a stray base **under the
corrected condition** (test 3 — the product comparison, stated so it survives the shares being
untied). *Depends:* F1, E3. *Source:* spec §4.2, §4.3; arch §4.1.

**F3. `ClassicEmissionOracle` + the cross-model check.**  ☐
Model B ported **test-only**, exactly as production keeps it
([`classic.rs`](../../../../src/ssr/cohort/read_model/classic.rs)) — an independent
implementation worth more as an oracle than as an alternative. Test: A and B agree on genotype
*ordering* over a fixture grid of (observation, candidate) pairs at matched parameters — the
independent-implementation check behind the whole path. *Depends:* F2. *Source:* spec §9; arch
§4.1.

> **Checkpoint F:** the seam holds two models, the default proven against the independent one.
> Pause for review.

### Milestone G — the censored term

**G1. `censored_emission`.**  ☐
Spec §5.2: the factorised form on pure candidates — the letter match on the witnessed prefix
times the closed-form tail `P(length ≥ ℓ | a)`, both geometric tails capped at E2's cutoffs —
and the exact sum over reachable stretchings on interrupted candidates. Tests: **the complement
identity** — `P(≥ ℓ) + P(< ℓ)` equals the truncated distribution's own total (one minus E2's
reported loss), **not 1** (spec §12 test 12); where the constraint admits exactly one length
change, censored equals complete **bit for bit**; **a partial never out-discriminates a
complete observation** on a stated parameter set (test 13). *Depends:* F2, E2. *Source:* spec
§5.2; arch §4.1.

### Milestone H — the STR row

**H1. The row, cache, and two-term mixture.**  ☐
`ssr.rs`: `genotype_log_likelihood_row<Model>` — emissions cached per
`(observation, candidate)` and reused across genotypes (**the cost is
`observations × candidates`, not `× genotypes`; that is the design, not an optimisation**);
complete and partial observations routed to `emission` / `censored_emission` by `ReadWitness`
through the guard iterator; the copy-weighted mixture with the outlier term,
`DEFAULT_OUTLIER_WEIGHT = 0.01` inherited and declared inherited. Tests: the junk term cancels
bit-for-bit for a read nothing explains (spec §12 test 6); ploidy generality with the corrected
between-ness claims (test 7 — pin both the matching case and the split case); order independence
(test 8); an instrumented count pinning one row build at `observations × candidates` emission
calls. *Depends:* F2, G1. *Source:* spec §2.1, §4.5, §8; arch §4.1.

**H2. The per-locus outlier spread + STR contamination.**  ☐
`reachable_length_count` computed from the candidate set and the cutoffs alone — **no cohort in
it**, the decided repair of production's cohort-wide `D`
([`em.rs:393`](../../../../src/ssr/cohort/em.rs)); and the three-term form
`(1 − λ − c)·copy-mixture + λ·uniform + c·seed(o)` with `SsrContamination` (`None` ⇒ the
two-term form), the seed's length distribution handed in from the prior's
`seed_length_distribution` ([`calling_prior.md`](calling_prior.md) E2). Tests: **the per-sample
property** — adding an unrelated sample's observations to the locus's cohort changes this
sample's row by zero bits; the `c = 0` ulp identity on this path too (the same named tolerance
as C1). **Own commit, do not bundle** — the cohort-dependent floor was a silent tenfold
difference between one sample and a panel; the zero-bit test is the oracle. *Depends:* H1;
prior plan E2. *Source:* spec §4.5, §4.5.1; arch §4.1.

> **Checkpoint H:** both rows are complete, pure, and proven; the loop plan can consume them.
> **Step 7 is complete as a set of pure functions.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | type-level compilation as the sibling plans' import surface; contract documented |
| B | **bitwise aggregation identity**; **production differential** — `standard_log_likelihood` reconciled term-by-term through the two recorded changes |
| C | the ulp-tolerance `c = 0` identity (named constant); a hand-computed contaminated case |
| D | hand-built compatibility fixtures — no-information, error-charge, pooled-verdict |
| E | existing stutter tests green through the rename; the reported-loss identity; the sums-to-one tripwire |
| F | the three ported property tests; **`ClassicEmissionOracle`** cross-model agreement |
| G | the complement-vs-truncated-total identity; the single-length bitwise case; the never-out-discriminates bound |
| H | junk-cancellation bit-for-bit; ploidy properties; the emission-call count; **the zero-bit cohort-independence test** |

## Out of scope (next plans)

- **Building contexts per locus from `StratumFits`, the `Lg` table, and row reuse across
  passes** — [`calling_loop.md`](calling_loop.md).
- **The change measurements** (spec §12 items 14–19: the coefficient, both end-to-end
  regressions, STR contamination on/off, zero-slippage path comparison, per-locus re-fit) — need
  genotypes; recorded in [`calling_bakeoffs.md`](calling_bakeoffs.md)'s out-of-scope with this
  plan named as the producer of the seams they measure.
- **`q_sum_other`'s producer** — candidate selection's, once that spec exists (see Scope).
- **The deferred items with homes** — per-read mismapping mixture, overdispersion, per-base STR
  rate, part-repeat estimator, purity adjustment, longer-than-read alleles (spec §10).
