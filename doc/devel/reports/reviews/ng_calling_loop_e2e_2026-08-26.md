# ng calling loop — E2e review: what three agents found

**Step:** E2e of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the repeat tract's prior
seed, read off the fit rather than constructed.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.
**What was reviewed:** the uncommitted working tree at `424a0808` plus the step's diff, in three
detached worktrees created from that commit, one brief each — arithmetic and control flow; tests and
mutation; design conformance and claim-checking.

---

## 1. What the three found, in one paragraph each

**Arithmetic.** No defect in the seed's own arithmetic: every edge the brief named — a candidate
below the reference by more than the span, a repeat count of zero, a three-class spectrum, a
one-candidate locus, the interaction of the two floors — was checked by running, and none produces a
`NaN`, a negative entry, an index panic or an entry below the concentration floor. **What it found
instead was a door.** `LengthSpectrum`'s fields were public, so the checks at the gather were not
the only way in: built by hand with an empty `weights`, `allele_span()` computed `(0 − 1) / 2` —
a subtract-with-overflow panic in debug and, **in release, a span of −1**, against which no candidate
is ever in reach, so every tract of that stratum came back with a degenerate `[1e-12, …]` prior and
nothing downstream could see it. An even class count did the same more quietly. It also found that
`SsrFitConfig::allele_span = 0` is a legal configuration the fit accepts and the gather then aborts
on, with a message naming a class count rather than the knob.

**Tests and mutation.** Five surviving mutations that mattered, of which **the first is the one to
learn from: the middle rung — the whole new fitting cost of this step — was never seeded from its
own weights.** Routing `PeriodsPooledTracts` to the flat arm left all 55 tests green, because the
one fixture that ever passed a pooled spectrum asserted only *share equals seed entry over seed
total*, an identity that holds for any shape. It also found the even-class fixture `[0.5, 0.5]`
failing *both* halves of its check, so deleting the odd-count condition changed nothing; a
`close()` helper whose tolerance never fell below `1e-12` — the size of the very floor it was used
to check — so the concentration-floor test passed with the floor deleted, in release *and* in
debug; a longer-than-the-candidate-set buffer claimed to be refused and never tested; and
`FrozenParameters::ssr_length_spectrum_at` with no test and no caller anywhere in the tree.

**Design conformance and claim-checking.** **12 wrong claims of about 95 checked, and every counted
figure was right again** — the failures were all statements about a mechanism or a location. The
Blocker was a supersession that only one side announced: `spec/calling_priors.md` §5 and
`arch/calling_priors.md` §5 still presented the deleted interface as *the* interface, with no banner,
while `population_diversity.md`'s own opening sentence — *"nothing has been built against §5's
version"* — was false, because E1 had built it and E2e had just deleted it. It also produced the
full list of retired sentences across seven documents, and independently reproduced the step's
measured cost ratio on its own machine (0.574 / 0.653 / 0.627 against the report's 0.60).

## 2. The findings that changed code

| finding | what it was | what was done |
|---|---|---|
| **Blocker** — the supersession was one-sided | two design documents still specified the deleted `fill_ssr_seed`, and the superseding document's own opening claim was false | banners on both §5s and on the file headers; the false sentence replaced **and the correction recorded**, because the half that was wrong is what a reader would have used to judge how much code was at stake |
| **Major** — `LengthSpectrum`'s fields were public | an empty or even-length `weights` reached the seed builder and produced a degenerate prior, silently in release | the type is a **struct with private fields** in a nested module, behind `fitted` and `stated_flat`; both check the class count and the concentration. *(It reads as three variants and is not an enum: an enum's variant fields carry the enum's visibility, so three variants would have left every check optional.)* |
| **Major** — the middle rung was never seeded from its own weights | routing the pooled rung to the flat arm left all 55 tests green | `the_pooled_rung_is_seeded_from_its_own_weights`, a twin of the top rung's test that pins four numbers |
| **Major** — the even-class fixture failed both halves | `[0.5, 0.5]` is short *and* even, so the odd-count check was untested | the fixture is `vec![0.25; 4]`, which only the odd check refuses |
| **Major** — `close()` had an absolute floor the size of the floor it checked | the concentration-floor test accepted zero | `close` is relative with no floor, and the floor is asserted with `assert_eq!` — it is exactly representable |
| **Major** — a longer buffer was never tested | both mis-sized-buffer tests passed a *shorter* one | `a_buffer_longer_than_the_candidate_set_is_refused` |
| **Major** — `ssr_length_spectrum_at` had no test | `+ 1` on its key passed all 55 tests | `the_run_answers_a_tracts_prior_shape_from_its_own_stratum`, on two strata differing in period, shape and concentration |
| **Major** — `SsrFitConfig::allele_span = 0` | a legal configuration that aborts the gather three modules later | refused in `fit_pooled`, where the message can name the knob, with a test |
| **Minor** — the pool's period was carried twice and unchecked | a pool filed under another period's key seeded that period's tracts from this one's spread | `assert_eq!` and a test |
| **Minor** — `strata_pooled` counted strata with no tracts | it says "contributed tracts to it" | counts members with a spanning read |
| **Minor** — the pool was cloned before the floor was applied | a period about to be discarded still paid for a full copy of its evidence | the floor is a sum over members, tested first |
| **Minor** — the floor's unit was untested | every fixture drew reads at every tract, so `tracts_with_reads()` and `tracts.len()` were one list | a fixture with four unread tracts |
| **Minor** — the pooled fit's `homozygote_excess` was unexercised | passing zeros changed nothing any fixture asserts | the same tracts read as a selfing panel and as an outbred one, **4.9% apart** in concentration |
| **Minor** — `converged`, the counters and the normalisation tolerance were unasserted | `.len().min(1)`, `converged: true` and a tolerance of `1e-2` all survived | four tests |
| **Minor** — `SHAPE_FLOOR` was only ever compared against itself | every assertion was written `12.0 * SHAPE_FLOOR` | a value pin, against production's `G0_FLOOR` |
| **Minor** — `should_panic` strings that did not discriminate | the pooled one matched all four of its function's messages | both spectrum-refusal tests now match the specific sentence, and the stratum-side one is at period 3, 17 repeats so its two interpolations cannot be swapped |
| **self-found** — `LocusInference::seed_diversity_unreachable` | a public output field whose whole subject was the deleted refusal | replaced by `length_spectrum_rung`, the carrier §8's third check needs |

## 3. The findings that changed only prose

Eleven, all from the claim-checking brief and all recorded in the implementation report or at the
site: a `tests::` back-reference naming a test that does not exist; "the three sit within a screen
of each other in a tract's parameter assembly", where the third has no caller at all; a goal claimed
met that §4 of the same report admits is not; "one pass over the candidate lengths" where there are
three; a scratch harness cited by a path not in the tree; a citation of §4.4 for a goal that lives
in §1; `ALPHA_REF` called a *share* two paragraphs after the module insists a share and a count of
chromosomes are different quantities; and the shipped doc comment that said the pooled fit *"roughly
doubles"* the repeat-tract half of a run, where the report's own measurement of the same thing —
made by the same author, in the same hour — says 60%.

**And one whole class of them:** sentences across seven documents that still describe the retired
construction. `spec/` and `arch/calling_priors.md`, `impl_plan/calling_prior.md`, all three
`candidate_alleles_ssr` documents, `spec/read_likelihoods.md`, `arch/read_likelihoods.md`, and one
test fixture's doc comment in `likelihood/ssr.rs`. Each now carries a banner or a corrected sentence.

## 4. What was found and deliberately not fixed

**`population_diversity.md` §5's *absent-or-present as a whole*.** The spec decides the tract
parameters are carried absent or present as a whole and that a tract in a run with no repeat-tract
parameters is refused by name; the implementation makes absence unrepresentable — `FrozenParameters`
holds a bare `&StratumFits` — and `a_run_that_fitted_no_stratum_states_the_constant` pins the ladder
answering instead. That is §4.4's *always answers* against §5's *refuse*, and the two are not
reconciled. **Left as a stated gap** rather than resolved here: giving `FrozenParameters` an optional
bundle changes a shared type and reaches beyond this step, and the consumer that would do the
refusing is E3b's tract branch.

**A candidate set every one of whose lengths is outside the fit's span** gets a prior with every
entry at the floor — a near-degenerate Dirichlet rather than the flat one the floor's documentation
implies — and the sibling export disagrees with it there, returning `1/K`. Unreachable in practice
because the reference candidate sits at offset 0, but nothing in `fill_ssr_seed` requires that, and
candidate selection can cut the reference allele. **Banked for E3b**, which is where the invariant
would be stated.

**`RepeatCount` carries both meanings.** `ssr_length_spectrum_at(period, reference_repeats)` and
`ssr_substitution_rate_at(read_group, period, candidate_repeats)` take the same type with opposite
required meanings, and the only guard is the parameter name — whose cost the step measures at 0.595
of the prior's mass moving to 0.091. A `TractRepeatCount` newtype would close it. **Not taken here**:
it is a shared-vocabulary change with callers on the other side of the seam.

## 5. What the reviews cost and what they were worth

Three agents, three worktrees from `424a0808` with the working tree applied as a patch, deleted
afterwards. **Twenty findings changed code and eleven changed prose.** Of the code findings, six
were tests that could not fail — which is the sixth consecutive step of this plan on which the
review's largest category was that, and the seventh time the accident was *every fixture shares one
coincidence*. This step's coincidences were: only one fixture ever used the middle rung, and its
assertion was an identity; the even-count fixture was also too short; and a floating-point
comparison helper's tolerance was the size of the constant it was used to check.
