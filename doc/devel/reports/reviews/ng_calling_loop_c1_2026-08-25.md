# Code Review: ng calling loop — C1, the flat first pass

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step C1
**Implementation report:** [ng_calling_loop_c1_2026-08-25.md](../implementations/ng_calling_loop_c1_2026-08-25.md)
**Fixes applied:** [fixes_applied_2026-08-25_v5.md](fixes_applied_2026-08-25_v5.md)

## 1. Scope

C1's working-tree diff on top of `abaebb12`: `PassPrior`, the `match` in `score_one_sample`, the
test helper's `flat_first` and explicit `starting_cohort`, and three new tests. One file,
+447/−92. **Two agents**, each in its own worktree — reliability + errors, and naming +
idiomatic + smells + refactor_safety + module_structure + step 8a. Two rather than three, in
proportion to a one-file diff; recorded rather than silent.

## 2. Verdict

**Request changes** — 0 Blockers, 8 Majors, 5 Minors. One Major is a defect in the function; two
are wrong claims in the diff's own prose; one is a design gap the spec does not answer.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --all -- --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | no warnings |
| `cargo test --lib` | 0 | `4599 passed` at review time (B2 left 4,596) |
| `cargo test --release --lib ng::calling --all-features` | 0 | green, as the CI step requires |

## 4. The defect

**The flat arm bypassed two release-held shape checks the seeded arm got for free.**
`prior_row`'s width was `PriorRow::new`'s to check and `sample_expected_copies`' was
`fill_sample_concentration`'s — and a flat pass enters neither. Measured on a 3-genotype locus
given a 2-entry prior row: the seeded arm panics in release, and the flat arm returns
`posterior = [0.199, 0.399, 0.402]` against the right `[0.25, 0.5, 0.25]` **in silence**, the
tail entry being the stale value the buffer arrived holding, carried through the normalisation
because the score loop's `zip` stops at the shortest row while the normalising loop divides all
of them. At a 3-allele locus with a 2-entry copies row the flat arm returned `[0.667, 0.667]` —
a diploid sample carrying 1.33 copies of a genome.

This is the module's declared failure mode reversed on one arm. `SampleScoringBuffers`' own doc
states the invariant that had become conditional.

## 5. Two claims of mine that were wrong

**The window was reported as though cohort size were not an axis.** The sweep covered advantage
× samples (3, 6, 20, 63) and its own output shows the two starts diverging at **0.5 nats for 20
and 63 samples**. The comment said there was "nothing to lose" at 0.5 nats and that the effect
bit near 1 nat "and nowhere else". That is `CLAUDE.md`'s named failure — a figure measured in
one corner, reported as a property of the caller.

**The rejected design's failure mode was wrong in both halves.** Handing the seeded path a first
pass's scratch does not panic on the cohort check and does not fall through to the bare seed in
release: it panics on *this function's own* release-held check that the sample's own copies are
finite — **in release as well as debug**. The bare-seed fall-through is real only when the
sample's own row is finite and the cohort row alone is `NaN`; probed separately, that gives
`concentration = [1.0, 0.5]` against a seed of `[1.0, 0.5]`.

## 6. The trap is a delay, not a different answer — and spec §3 says otherwise

The largest finding, and it corrects the step's own justification.

**Measured at 63 samples**, alternative copies per sample: the seeded start sits at 0.151 at six
passes, **0.633 at nine — where it flips to heterozygous** — and both starts reach 0.767332 by
thirty. At three samples the flip is between ten and sixteen passes; both reach 0.549452.

**So both starts reach the same fixed point.** What the flat pass buys on this fixture is
**eight passes**, not a different call. Spec §3 says the seeded loop *"converges, and it
converges to no-variant, having never let the reads speak"* — **that stronger claim was not
reproduced.** A rare-variant shape was tried too, a handful of carriers among 60 samples whose
reads are firmly homozygous reference, swept over carrier count (1, 3, 6) and advantage (1, 2, 4
nats) at 50 passes: the two starts agreed in **every** cell.

**This does not undermine C1.** The mechanical reason for the flat pass — the cohort's copies do
not exist on the first pass — is unassailable, and the choice is inherited from GATK and from
production. What it changes is what the code may claim: the delay is real and matters against a
cap of 50 passes and production's observed 3-to-5-pass convergence, but *permanent* no-variant
convergence is not something this work could demonstrate. **Raised for the owner against §3.**

## 7. The design gap: a silent sample votes for a 50% allele frequency

Spec §7 says a sample with no coverage *"scores every genotype alike, so the prior decides it
alone — the right answer rather than a special case"*. **On a flat pass there is no prior to
decide it.** Its posterior is the normalised genotype table and its expected copies come out as
the average genotype — a full copy of the alternative at a biallelic locus, the same
contribution as a confident heterozygote, and about 1,000× what the seeded start gives it.

Measured at 63 samples: three silent samples put **3.0027** alternative copies into pass 1
against 0.0030 with none; six silent samples still leave 0.48 after two passes against 0.000001.
At three reads a position roughly one sample in twenty is silent at any given position, so this
is the ordinary case rather than a corner.

**Neither §3 nor §7 says what should happen.** Pinned by a test rather than changed, and raised
for the owner.

## 8. Answers to the brief's questions

- **`reads_the_cohort` is consulted in exactly one place and is correct** — no seeded pass
  escapes the own-copies check; mutations forcing the predicate either way are both killed.
- **The zero prior row is equivalent to no prior** for every degenerate likelihood row that could
  be constructed. A non-zero constant also survives, differing by one ulp, which the prior seam's
  contract explicitly permits.
- **Neither stale buffer can reach an answer**: `prepare_for_locus` refills both per locus,
  `fill_sample_concentration` rewrites the concentration whole, and both shipped priors write the
  per-allele workspace before reading it.
- **At one sample spec §7's claim holds bitwise** — flat pass 1 plus seeded pass 2 gives winner
  `[0]` and cohort `[1.998642703842511, 0.0013572961574892799]`, identical to a single seeded
  pass. **At two samples the flat start does not hold the heterozygote either**, so "the flat
  start finds the variant" is a property of three samples and up on this fixture.
- **`PassPrior` belongs where it is** — one consumer, file-local, and the plan puts D1's loop in
  the same file. It is a per-*sample* value under a per-*pass* name, though: writing D1's three
  nested loops shows the per-sample form compiles and the hoist the name invites fails with
  `E0425`.

## 9. Step 8a — 15 claims, 11 correct

The ones easiest to get wrong were right: the prior row prints back as
`[0.6931471805599453, -6.907755278982137, -7.600402584500431]`; the gap is `7.600902`, and
`exp(7.600902) = 1999.9999999999998`, so "about 2,000 to 1" is exact; 1 nat is `4.3429` Phred;
`(1 − 3θ/2)/θ` at `θ = 0.001` is `29.993` Phred; one flat pass leaves `0.7310585786300917` copies
per sample, which is exactly `1/(1 + e⁻¹)`; twelve copies over six samples came out as `12` to
the bit. The three wrong ones are §5 and §6 above.

## 10. Mechanical

**`CohortFixture` was inserted between `run_passes`'s doc comment and `run_passes`**, so all 29
lines of that doc documented the struct and `run_passes` had none. Nothing warns — both are
private test items.

## 11. What's good

The shape of the decision. Making the pass's prior a value rather than a code path is what
`arch/calling_em_loop.md` §2.1 asks for, and both agents independently concluded `PassPrior`
belongs in the file it is in. And the counterfactual mutation — replacing the flat arm with the
seeded one — fails exactly the three tests written for it and nothing else, which is the cleanest
oracle any step in this milestone has had.
