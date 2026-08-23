# ng genotype prior — E1: the repeat-tract seed, and the total that reproduces the measurement

*Implementation report, 2026-08-23. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step E1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone E, on top of `0b019e0d`.*

## 1. What it is

`fill_ssr_seed` writes the concentration the genotype prior starts a repeat tract from — one
positive number per candidate length, read as *chromosomes the prior behaves as though it had
already seen*. Two independent parameters set it, and separating them is the one place ng departs
from production's shape rather than porting it
([`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §5.1):

```text
where the mass sits   w_j ∝ max(decay ^ |repeat count of j − cohort modal count|, 1e-12),  Σ w = 1
how much there is     A = D / (1 − c − D),   c = Σ_j w_j²    (the shape's Simpson index)
the seed              α_j = A · w_j,  floored at MIN_ALT_CONCENTRATION
```

The **shape** is production's `G₀` ported unchanged
([`allele_freq_prior.rs`](../../../../src/ssr/cohort/allele_freq_prior.rs)). The **total** is new,
and it is the whole point of the step.

## 2. The units error the total exists to avoid

Gene diversity `D` is a probability — the chance two copies drawn at random carry different
lengths. A concentration is a count of chromosomes. Setting the total to `D` equates the two, and
spec §5.1 records that the design document said exactly that until 2026-08-19. What a Dirichlet
with total `A` and Simpson index `c` implies is `A(1 − c)/(A + 1)`, so a total of `D` asserts
`D(1 − c)/(D + 1)`, always less. On 1,236 polymorphic tomato repeat tracts at the fallback decay
the prior would assert a median 0.40 of what was measured, tenth percentile 0.22.

Inverting the identity gives `A = D/(1 − c − D)`, which is what the code uses, and
`the_seed_implies_the_diversity_that_was_measured` is the oracle: it reads `A` and `c` **back off
the buffer that was written** and checks `A(1 − c)/(A + 1)` against the measurement, so it cannot
agree with the builder merely by re-running the builder's own algebra. Worst relative error over
the sweep's 118 seeded cases: **1.12e-15**, about seven units in the last place of a diversity of
0.087.

One reviewer derived the identity independently two ways — the Pólya urn, and `E[p_j²]` — and both
give the code's expression, so the mathematics is confirmed rather than merely restated.

## 3. The refusal, and the range it has to work across

`A(1 − c)/(A + 1)` rises to `1 − c` and stops, so a tract whose measured diversity is at or above
`1 − c` has no total that reproduces it. `SsrSeedOutcome::DiversityUnreachable` reports it with
both numbers and never rescales silently (spec §12 test 11, and open as spec Q2).

**How often it fires was the review's biggest finding, and it is a fact about the panel rather
than about the caller.** The one-in-ten figure the spec quotes — 119 of 1,236 tomato tracts at the
fallback decay — is a 63-accession selfing crop at about three reads a position. At the other end
of the committed cohort range the refusal is the rule:

| panel | candidate lengths a tract shows | ceiling the shape can imply | measured `D` | refused |
|---|---|---|---|---|
| tomato, 63 accessions | several | varies | median 0.087 | 119 of 1,236 |
| one outbred genome | at most 3 | 0.444 (two) to 0.625 (three) | ≈ 0.72 (GIAB HG002) | every tract |

The pre-pass fits this quantity at every cohort size down to one
([`parameter_prepass_cohort.md`](../../ng/spec/parameter_prepass_cohort.md) §3), so at one sample
it returns that genome's own repeat diversity — about 0.72 on HG002, where 72 tandem repeats in
100 are heterozygous (spec §5.3). No decay rescues it: the ceiling saturates at `1 − 5/27 = 0.815`
at the fallback decay however many lengths a tract carries. Pinned by
`a_single_outbred_genome_is_refused_at_every_tract`, and written into spec §5.1, §12 test 11 and
Q2, because **whichever policy Q2 settles has to work at ten refusals in ten**, which rules out
any candidate that is affordable only because it is rare.

**One locus the rule does not reach.** A tract with a single candidate length has a Simpson index
of exactly 1 and therefore a ceiling of 0, so the rule would refuse every monomorphic tract at any
measurement including zero. There is nothing to refuse — one length is one genotype, whose prior
probability is 1 at any positive concentration — so it is seeded at `ALPHA_REF` and the rule
starts at two lengths. Recorded in arch §5 and spec §5.1, which both stated the rule without it.

## 4. Departures from arch §5, all three recorded there

- **Named `fill_ssr_seed`**, matching the module's other buffer-fillers.
- **The two scalars are checked types in `ng/types.rs`** — `RepeatGeneDiversity` and
  `SeedDecayPerRepeat` — with `try_new` returning `DomainError`, beside `ExpectedHeterozygosity`.
  Arch sketched bare `f64`. They are the pre-pass's outputs like every other measured scalar the
  caller consumes, so a degenerate fit returning a `NaN` is a run to refuse with a message rather
  than a process to abort. That is why they are **not** in `genotype_prior`'s `checked` module,
  whose constructors panic. The fallback decay is `SeedDecayPerRepeat::FALLBACK`, an associated
  constant following `ExpectedHeterozygosity::SPECIES_FALLBACK`, rather than a free-standing
  `DEFAULT_G0_FALLBACK_DECAY` — as a loose `f64` it is exactly as constructible into a stutter
  one-step share as into this, which is the trap the rename was for.
- **The refusal withholds the concentration** and hands back the buffer instead, where arch and
  the plan both describe a marker on a returned value ("the loop uses the ceiling total"). The
  difference from `SpectrumMatch` on the SNP/indel path is deliberate and now written down in
  arch §5: the spectrum fit runs **once per run** and the run cannot start without a seed, so
  withholding would leave the caller nothing to do but invent one; the STR refusal is **per
  locus** and the loop can pick a policy and carry on.

## 5. What the reviews found

Four reviewers, each in its own worktree, each given the gate output rather than re-running it.

**Three should-fixes on this step, all applied:**

- **The shape floor was live code no test pinned.** Two reviewers found independently that
  deleting `.max(SHAPE_FLOOR)` left every test green. The test named for it asserted only that a
  far candidate's *concentration* was positive — which `MIN_ALT_CONCENTRATION` guarantees on the
  seeded path whatever the shape did. It now asserts the far candidate's **share**, exactly, and a
  second test covers the case the floor's own doc names: a tract whose every candidate sits past
  the underflow distance, which without the floor divides by a total of zero and fills the buffer
  with `NaN`.
- **The outcome type's central safety claim was false.** It said the concentration was reachable
  only through `Seeded`; `Concentration::new` is public, and a caller can wrap the shape by hand.
  Closing that would break the design, because arch §5's provisional policy for these loci needs
  exactly that constructor. The claim is replaced by what the type actually delivers: the mistake
  stays representable, but it becomes an explicit decision with a name on it rather than something
  a caller falls into by reading a buffer it was never handed.
- **`fill_seed_shape` had no length check**, safe only because its one caller asserted first — and
  plan step E2 makes it public. The check moved into it.

**Six wrong numbers in this step's own doc comments, every one a claim about its own fixture** —
the same failure the three preceding steps produced. Corrected: "five units in the last place" was
seven; "3 in a trillion" was 1.8; "three passes" was four; "three lengths checked" was two; "more
than ten thousand chromosomes" was about a million; and the underflow offset was given as both
"about 1,075" and "past 1,074" in the same commit. A seventh claim — that the seven refusals in
the identity sweep were "the steepest decays against the widest diversities" — was wrong about the
mechanism: all seven sit at one diversity, and five are the two-length tract at **every** decay
including the flattest, because two lengths cannot imply more than 0.5 however flat the shape is.

**Four claims about other files were wrong and are corrected:** the shape floor's reason was
attributed to production, which gives only the underflow half; the two floors were said to
coincide because both are "any representable positive number", which is true of production's and
not of `MIN_ALT_CONCENTRATION`, whose value is sized; and the comparison with the stutter model's
one-step share was wrong in three ways — the stutter geometrics are success probabilities rather
than decays, the share is of reads that slipped by whole repeats rather than of all copying
errors, and they are clamped to `[0.01, 0.99]` rather than `(0, 1]`.

**Mutation testing, four reviewers.** Twenty-nine distinct mutations across the step and the
commit under it. Survivors on this step: the shape floor (fixed above, now caught), `#[must_use]`
on the outcome (a compile-time warning no test can observe — recorded, not fixed), and the
saturating integer clamp on the offset, whose wrapping alternative needs two repeat counts 2.1
billion units apart and so cannot be reached through the function's arguments. Two mutations were
proved to be no-ops rather than coverage gaps: `headroom <= 0.0` against `headroom < 0.0`, and
dropping the finiteness guard — both indistinguishable over 32,200 boundary inputs, for the reason
the code already gives.

## 6. What the tests pin

Nineteen, up from thirteen before the review.

| test | what it holds |
|---|---|
| `the_seed_implies_the_diversity_that_was_measured` | spec §12 test 10 — the identity, read back off the written buffer, over 118 of the sweep's 125 combinations, with the other 7 counted as refusals so a shrinking sweep cannot hide |
| `the_concentration_floor_lifts_the_implied_diversity_and_by_how_much` | the floored path is a separate measurement: upward, bounded, 3 parts in 10 million at 50 lengths and a diversity of 1 in 10,000 |
| `a_diversity_the_shape_cannot_hold_is_refused_exactly_at_the_bound` | spec §12 test 11 — four spreads at three decays, each bound found from the shape |
| `a_refusal_hands_back_the_shape_and_asserts_no_total` | the buffer comes back through the outcome, summing to one |
| `a_single_outbred_genome_is_refused_at_every_tract` | the cohort end of the range, and the 0.815 saturation |
| `a_candidate_too_far_for_the_decay_to_reach_keeps_its_share` | the shape floor, measured as a share so the concentration floor cannot stand in for it |
| `a_tract_whose_every_candidate_is_past_the_underflow_distance_is_still_seeded` | the floor's second job — the normalisation |
| `an_offset_past_the_integer_range_is_clamped_rather_than_wrapped` | wrapping would invert the shape onto the far candidate |
| `a_second_spelling_of_one_length_raises_the_ceiling` | spec Q3 with a size: 0.444 to 0.625 |
| `a_locus_with_one_candidate_length_is_seeded_rather_than_refused` | the exception, and that its ceiling really is zero |
| `a_cohort_with_no_repeat_variation_seeds_every_candidate_at_the_floor` | the `θ = 0` twin, and the 1.8-in-a-trillion crossover where the geometry returns |
| `the_total_climbs_without_bound_as_the_measurement_nears_the_ceiling` | the pole, and that the seed reaches a million chromosomes at a gap of one in a million |
| `the_largest_total_any_shape_can_ask_for_is_finite` | `2^53 − 1`, so the overflow guard is unreachable by construction |
| `mass_falls_off_by_the_decay_...`, `two_candidates_of_one_length_...` | the port's fidelity |
| the four refusal tests | both length checks in both fillers, and the two scalars' domains |

## 7. Open items this step touched

- **`arch/calling_priors.md` §5 cited `read_likelihoods.md §8` for the provenance channel both
  refusal markers travel on. That section does not exist** — the arch sibling ends at §7. Corrected
  to §1.4, which is where provenance actually lives. This was on the inherited open-items list and
  is now closed for this citation; the other four `§8` references in that file point at the **spec**
  sibling, whose §8 does exist.
- **Doc links to `#[cfg(test)]` test functions do not resolve**, which a reviewer found on the
  first draft of this step. They are written here as plain code spans instead, so this step adds
  none. `cargo doc --no-deps` on the branch now reports **12 warnings, all "redundant explicit link
  target", and no unresolved links at all** — so the handoff's note of "17 pre-existing unresolved
  intra-doc links" is stale and there is nothing outstanding to inherit.
- **Spec §5.1 says production "carries the question as a deferred re-tune"** with a quoted label
  that appears nowhere in production. Flagged for whoever owns the spec; the code quotes the spec
  faithfully.

## 8. Gates

Green in the container: `cargo fmt --check`, `cargo clippy --lib --tests --all-features -D
warnings`, `cargo test --lib`.
