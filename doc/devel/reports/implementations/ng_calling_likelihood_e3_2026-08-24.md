# ng read likelihoods — E3: seven shares from three fitted numbers, and the tripwire that catches three silent failures

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step E3 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone E, on
top of `43bbb149`. **This completes Milestone E — Checkpoint E.***

## 1. What it is

The parameters fit produces **three** numbers per read group per stratum
([`Slippage`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)): how often a read slips at
all, which way it slips when it does, and how fast bigger slips fall away. The stutter
distribution wants **seven**. `stutter_rates_for` is the conversion, and it is now the only place
in ng that performs it.

Alongside it, spec §12's fourth test: **the distribution sums to one over its whole support**.

## 2. The complement, which is the whole reason this function exists

`Slippage::fall_off` is the probability of **carrying on** to a further step — its own
`read_probabilities` weights step *n* as `(1 − fall_off) · fall_off^(n − 1)`. A `StutterRates`
one-step share is the probability of **stopping at one step**. **They are complements**, so the
conversion is `one_step_share = 1 − fall_off`.

Getting it backwards inverts the size distribution — large slips become the common ones — and
nothing crashes. It is the first trap spec §4.2 names, and production writes the same conversion in
the same words. Two tests fail if it is dropped, and the fixture uses a fall-off of **0.35
deliberately**: at a half the two are equal and the mistake is invisible.

## 3. How the level becomes four shares

| share | from |
|---|---|
| `whole_repeat_shorter_share` | `level × shorter_share` |
| `whole_repeat_longer_share` | `level × (1 − shorter_share)` |
| `part_repeat_shorter_share` | `level × 0.05 × shorter_share` |
| `part_repeat_longer_share` | `level × 0.05 × (1 − shorter_share)` |
| both one-step shares | `1 − fall_off` |

The part-repeat mass is **added on top** of the level rather than carved out of it, so the four
direction shares total `1.05 × level` — production's shape, reproduced deliberately.

**Two of the seven are placeholders and are named as such**, because spec §4.2 requires them
recorded rather than mistaken for estimates:

- **`PART_REPEAT_SHARE_OF_WHOLE = 0.05`** — production's `OUT_FRAME_REL`, whose own comment calls
  it "the Step-4 declared estimator … pinned to a real per-period estimate in Step 5". That
  estimate was never made. Every part-repeat score in ng is therefore a fixed twentieth of the
  whole-repeat one, per read group and stratum, whatever the motif.
- **The two one-step shares tied to one number**, which HipSTR keeps independent. §10 gives this
  one a home.

Both are pinned by tests, so replacing either is a visible change rather than a silent one.

## 4. The tripwire, and the three failures it catches

`the_distribution_sums_to_one_over_its_whole_support` sums `probability` over every length change
the model scores, for periods 2 to 6, direction splits from symmetric to five-to-one, two slippage
levels and three one-step shares — and requires one. **No production test pins this.** All three
of the silent failures spec §12 test 4 names were injected here and caught:

| the mutation | what it does | caught by |
|---|---|---|
| the one-step share read as its complement | the ten scored steps hold `1 − 0.95^10`, four tenths of a branch's mass, so six tenths of every slip vanishes | the tripwire, plus 1 other |
| the same-length share not the remainder (`+1e-4`) | the total moves by exactly that error | the tripwire, plus 7 others |
| the part-repeat geometric indexed by `Δ` rather than the compressed rank | its weights skip the multiples of the period and sum short | the tripwire, plus 6 others |

**Period 1 is excluded and has its own assertion**, because there the part-repeat branch is
unreachable by construction and the total is *supposed* to fall short by exactly that branch's
mass — which E2's `a_mononucleotide_candidate_loses_the_whole_part_repeat_mass` already pins.

A second test measures the complement's cost rather than asserting it: a complemented share leaves
the distribution short by **0.1257 at a slippage level of 0.2**, against the tripwire's tolerance of
1e-9 — eight orders of margin.

The tripwire also cross-checks E2: alongside "the total is one to 1e-9" it requires the total to
equal `1 − truncated_mass_lost(...)` to 1e-12. Two pieces of machinery derived completely
differently, agreeing.

## 5. One deviation: where `stutter_rates_for` lives

The plan puts Milestone E's three changes in
[`alignment/stutter.rs`](../../../../src/ng/alignment/stutter.rs). **The conversion is in
[`calling/likelihood/ssr_emission.rs`](../../../../src/ng/calling/likelihood/ssr_emission.rs)
instead** (a new file, which the plan's scope line already names as this plan's).

The reason is a module edge. `alignment` and `parameter_estimation` are siblings and **neither
imports the other** — checked, not assumed. `stutter_rates_for` takes a `Slippage`, so putting it
in the alignment module would make the shared *distribution* depend on the *fitting* — the
direction `StutterModel`'s own contract disclaims in the sentence "this type fits nothing", added
at E1. `ng::calling::likelihood` already depends on both, and spec §7 puts the reading of frozen
parameters on that side of the boundary.

The sums-to-one tripwire stayed in `stutter.rs`, where the distribution is.

## 6. Validation

| command | result |
|---|---|
| `./scripts/dev.sh cargo test` — library target | **4,376 passed, 0 failed, 14 ignored** (4,369 at E2) |
| `./scripts/dev.sh cargo test --all-features` — every target | 4,470 / 0 / 18 |
| `clippy --lib --all-features --tests -- -D warnings` | exit 0 |
| `cargo fmt --check` | exit 0 |

**26 tests** in `ng::alignment::stutter` (24 at E2) and **5** in
`ng::calling::likelihood::ssr_emission`. Each of the three mutations in §4 was injected and
restored from a byte-compared copy.

## 7. Deferred for the owner

Carried forward from E2, plus one new:

1. **Spec §4.2's prose and its figures disagree** about whether a tract can lose its last repeat
   (E2's report §3). This branch follows the figures.
2. **`stutter_rates_for` has no production caller yet** — F1's `SsrScoringContext` is the consumer.
   It is exercised only by its own tests.
3. **The plan says Milestone E's three changes are in `alignment/stutter.rs`**; one of them is not,
   for the module-edge reason in §5. If that reasoning is wrong, moving it is a two-line change.
