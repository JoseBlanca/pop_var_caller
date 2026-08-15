# B1 — the depth ladder gains ten rungs

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md) step B1.
**Design authority:** [spec/parameter_prepass_generic.md](../../ng/spec/parameter_prepass_generic.md) §4;
[spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md) §2.2 and §4.1.
**Date:** 2026-08-14.

---

## 1. What changed

`DepthBinEdges::for_census()` — thirty bins where `new()` has twenty, topping out at **1,498** reads a
position where `new()` tops out at 124.

**It is one ladder at two lengths, not two ladders.** `for_census` takes the twenty-bin ladder's tops
as they are and appends ten more at the ratio those were built with, so every rung of the shorter
ladder is a rung of the longer one **by construction**. Rebuilding a thirty-rung geometric ladder over
the wider range instead would move the shared rungs by a base or two and quietly break that. The
ratio itself is now named once, `widening_ratio`, because two derivations that drifted in the last
decimal would round a shared rung differently.

The ten new tops: **159, 204, 262, 336, 431, 553, 709, 910, 1168, 1498.**

**Why the reach and not the resolution.** A position a sample carries two copies of holds about twice
that sample's median depth, so the signal dies wherever twice the median runs past what the encoding
can say. At a top of 124 that is 76 reads a position for a doubled position to stop reading deeper
than an ordinary one, and 98 for the two to be written identically — inside the range this caller
commits to (records spec §4.1).

**It costs nothing.** The stored code is five bits, holding 32 values; twenty bins plus the
never-walked sentinel used 21 and left eleven spare. Thirty plus the sentinel is 31. A compile-time
assertion in `census.rs` now ties the two together, so a ladder grown past what the field holds fails
to build rather than writing its top rung as the code for a bug.

## 2. The histogram route is untouched

`DepthBinEdges::new()` is unchanged and its cell table is still **583**, asserted in the same test as
the refinement. The plan's verification row asks for exactly that, and it is what "a census code maps
to a per-sample-route bin by collapsing everything above 124" means: the other route keeps twenty
rungs.

## 3. Tests

| test | what it pins |
|---|---|
| `every_edge_of_the_twenty_bin_ladder_is_an_edge_of_the_census_one` | the refinement, the ten new tops as literals, and the histogram route's 583 cells |
| `the_two_ladders_answer_alike_at_every_depth_up_to_the_shorter_ones_cap` | `bin_for` agrees on both ladders at every depth 0 to 124 |
| `a_doubled_depth_lands_at_least_one_rung_higher_all_the_way_to_the_ceiling` | **the property §4 leans on**, at every depth from 4 to half the ceiling — a ladder that saturated early would pass the first test and fail this one |

## 4. The census does not adopt the longer ladder here — it adopts it in B2

**A recorded deviation, and the reason is a coupling worth naming.** `recorded_depths` answers "the
cap and nothing else" for the ladder's **top bin**, because until now the ladder's top and the
per-position depth cap were the same number: a position deeper than 124 was thinned to exactly 124
before it was recorded, so the top bin held one depth rather than the twenty-seven it nominally spans.

Handing the census the thirty-rung ladder while the walk still clips at 124 breaks that. Positions at
124 would land in bin 19, which is no longer the top, so `recorded_depths` would return 98 to 124 —
a range where it used to return a point — and the fit's depth term would sum over twenty-seven depths
instead of one. On the trio at a few hundred reads a position that is most positions.

That is B2's subject exactly: the depth stops being clipped, and `recorded_depths` has to be revisited
with it. Adopting the ladder here would move numbers in B1 for a reason belonging to B2 and cost the
bisect the plan's second principle asks for. **So B1 is the ladder and its tests; B2 adopts it.**

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo test --lib ng::parameter_estimation::generic::depth_bins` | `18 passed; 0 failed` — 15 before, 3 added |
| `cargo test --lib ng::parameter_estimation::joint::census` | `24 passed; 0 failed` |
| the 88-second tomato oracle, against the A3 run | **byte-identical, 0 differing lines** — which is the point: nothing reads the longer ladder yet |
