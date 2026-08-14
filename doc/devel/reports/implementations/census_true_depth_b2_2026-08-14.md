# B2 — the depth recorded is the position's own

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md) step B2.
**Design authority:** [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§2.2 and §4.1; test §7.11.
**Date:** 2026-08-14.

---

## 1. What changed

**The walk records what the position actually held.** It used to record `min(depth, cap)`. The cap
still thins the allele counts beside it, proportionally, so the fractions they showed survive — that
is what bounds the sparse list at high depth and what keeps a count inside its own field.

**The census moved onto the thirty-rung ladder** at the same time, which is why B1 left it alone: the
two changes have one consequence between them and separating them would have moved numbers twice.

**`recorded_depths` is gone.** It answered `depth_range` except at the ladder's top rung, where it
answered the cap and nothing else. That was right only while the ladder's top and the cap were the
same number — and it was *wrong* even then for a position that genuinely held 100 reads, which was
never thinned and yet read back as 124.

**`DepthCap::denominator_for` replaces it, and it answers the other question.** Given the depths a
stored code stands for, it returns the depths the *counts* were taken out of — both ends clamped to
the cap. So a consumer asking *how many copies is this* reads the depth, and one dividing a count by
a depth puts the depth into the counts' units first. The cap comes from the recording terms, which
the fit has already refused to pool across.

Three consumers were pointed at it, each after asking which question it was really asking:

| site | wants |
|---|---|
| `EvidenceCursor::mean_depth` | the counts' denominator — it is the centre of the prior that weights the range the reference-read count is inferred over |
| `EvidenceCursor::advance`'s per-position range | the counts' denominator — the thinned counts are subtracted from it, and subtracting them from an unthinned depth charges reference reads the position never had |
| `contamination::markers` | the counts' denominator — it divides alternative reads by it |

**Nothing yet wants the true depth.** The consumer that will is the duplicated-copy class's depth
term (spec §4), and that class ships off. This step makes the number available; it does not spend it.

## 2. What moved, measured

Both oracles, before and after, at four container CPUs.

**Tomato — eight accessions from 2.4 to 30.6 reads a position: nothing moved at all.** The whole diff
is one header line in a diagnostic table. No fitted number, and not even the log-likelihood.

**The GIAB trio — three samples at a few hundred reads a position:**

| | before | after |
|---|---|---|
| read group 0, error rate at an ordinary position | 0.00491 | **0.00491** |
| …at a mismapped one | 0.0481 | **0.0481** |
| positions mismapped | 0.0083 | **0.0083** |
| carrying only the reference | 0.9983 | **0.9983** |
| only a non-reference base | 0.00020 | **0.00020** |
| the segregating density | Beta(3.459, 6.224) | **Beta(3.459, 6.224)** |
| expected heterozygosity | 0.638 /kb | **0.638 /kb** |
| HG002 heterozygosity | 0.786 /kb | **0.786 /kb** |
| HG003 heterozygosity | 0.601 /kb | **0.601 /kb** |
| homozygote excess, every sample | 0.000 | **0.000** |
| **log-likelihood** | −859,732 | **−858,973** |

**Every fitted parameter is unchanged to the precision the harness prints. The log-likelihood moved
by 759 nats.**

### 2.1 Why so little moved, and it sharpens the specification's prediction

The specification predicts change "confined to positions carrying more than 124 reads a position".
**Measured, the change is confined to positions carrying between 98 and 124** — a much narrower band —
and the reason is that above the cap the new code and the old agree exactly:

- **Above 124 reads** the walk used to clip the depth to 124, whose bin read back as `124..=124`. Now
  the walk keeps the true depth, whose bin `denominator_for` clamps to `124..=124`. Identical.
- **Between 98 and 124** the walk clipped nothing, but the old top-bin rule read every code in that
  bin back as exactly 124. Now it reads back as 98 to 124. **That is the correction.**

On the trio that band is about 1 position in 1,000: the occupancy table moves from 99.9% above the cap
to 99.8%, with the remainder crossing into the ranged column. On tomato there is no such position at
all, which is why nothing moved there.

So the step's visible effect on today's estimator is a small correction in a thin band. **Its purpose
is what it makes possible**, not what it moves: the depth a position actually held is now in the
census, where before it was destroyed at 124.

## 3. Tests

Two added, both spec §7.11, and both needed a fixture the module did not have — `writer_over_capped`,
because a cap *below* the ladder's ceiling is the only way to see the cap at all. Every earlier
writer fixture sat at 1 to 9 reads against a cap of 124.

- `a_position_above_the_cap_keeps_its_own_depth_and_thins_only_the_counts` — 40 reads at a cap of 20,
  a quarter of them non-reference. Asserts the code stands for 40 and not for 20, that the count came
  back 5, that `min(depth, cap)` is `20..=20`, and that the quarter survives.
- `a_position_below_the_cap_is_its_own_denominator` — the case every other fixture is in, so the
  clamp is shown not to be firing always.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo test --lib ng::parameter_estimation::joint::census` | `26 passed; 0 failed` (24 before) |
| `cargo test --lib ng::parameter_estimation::joint` | `71 passed; 0 failed` |
| `cargo test --lib` | `3575 passed; 0 failed; 11 ignored` |
| tomato oracle, 88 s | one header line differs; no fitted number |
| trio oracle, 74 s | §2's table |
