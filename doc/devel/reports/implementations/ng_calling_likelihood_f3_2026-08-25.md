# ng read likelihoods — F3: the second opinion, and where it disagrees

*Implementation report, 2026-08-25. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step F3 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). **This completes
Milestone F — Checkpoint F.***

## 1. What it is

`ClassicEmissionOracle` — Model B, ported **test-only**, as production keeps its own. It is worth
more as a second opinion than as a model anyone runs.

`Σ_n S(n) · avg_v align(observation | candidate ⊕ n repeats)`: the whole-repeat stutter mass
marginalised **outside** a sequence-versus-sequence align that absorbs whatever length is left
over as slop at the tract's ends.

## 2. What makes it independent, and what it shares

**Model A picks one length change** — the one the read actually shows — scores its probability,
and compares letters over sequences it has already made equal. **Model B sums over every
whole-repeat change** and lets the aligner explain the residual. Where A has an explicit
part-repeat branch with its own two fitted parameters, B has none.

So the two explain a read's length by genuinely different routes, which is what makes agreement
between them evidence. **What they share is the placement enumeration** — production shares it too
— so the independence is in how length is explained, not in where a slip may land. Stated rather
than implied, because an oracle that shared the thing under test would be worth nothing.

The align half is **composed, not written**: `SsrSequenceMarginal`
([`alignment/ssr_marginal_sequence.rs`](../../../../src/ng/alignment/ssr_marginal_sequence.rs)) is
already a faithful port of production's `align_subst`, end-gap slop and all.

## 3. What the check found

**Where the observation's length differs from the candidate's by a whole number of repeats — the
case Model B is built for — the two agree completely.** Across four observations against four
candidates the full ranking is identical every time, and the winning score agrees to better than
one part in ten thousand.

**Where it does not, they part company, and the divergence is the choice rather than a defect.** A
read whose length is not a whole number of repeats from the candidate is explained differently by
construction: A charges the fitted part-repeat share; B absorbs the odd base as slop at the flat
per-base rate. At ε = 0.001 the slop route is about **nine times cheaper** than the part-repeat
event, so the two prefer different candidates.

Measured, on a ten-base read against a three-repeat and a four-repeat candidate:

| | three-repeat candidate | four-repeat candidate | prefers |
|---|---|---|---|
| **Model A** | 1.094 × 10⁻⁴ | 1.869 × 10⁻⁴ | the four-repeat one |
| **Model B** | 9.706 × 10⁻⁴ | 1.167 × 10⁻⁵ | the three-repeat one |

**Model A is the one to trust here, and this is the reason it was chosen**: an explicit
part-repeat branch with its own fitted parameters is what the comparison behind spec §4.1 picked
it for. The test pins the disagreement with its sizes, so it stays a measured, deliberate thing
rather than a surprise the first time a ranking differs.

*(A first grid of pure candidates against exact-repeat observations agreed on everything, ratio
1.0000 at every winner. That grid could not have found this, and it is the reason the shipped one
carries part-repeat lengths, a substitution and an interrupted tract.)*

## 4. Validation

| command | result |
|---|---|
| `./scripts/dev.sh cargo test` — library target | **4,392 passed, 0 failed, 14 ignored** (4,390 at F2) |
| `./scripts/dev.sh cargo test --all-features` — every target | 4,486 / 0 / 18 |
| `clippy --lib --all-features --tests -- -D warnings` | exit 0 |
| `cargo fmt --check` | exit 0 |

**13 tests** in `ng::calling::likelihood::ssr_emission`, against 11 at F2.

*(One full-suite run failed first with a **rustc internal compiler error** in the codegen backend
while building `examples/ng_inbreeding_resolution`, which nothing here touches. It did not
reproduce — the immediate re-run was clean. This is the second time in this plan; the first was on
a different example. Recorded because a green log that follows a red one should say why.)*

## 5. Deferred for the owner

Unchanged from F2, and none of them blocks Milestone G:

1. **`unreachable_mass` understates the loss on interrupted tracts** — it sees a period and a total
   repeat count, where reachability is per run. The sharpest of these.
2. **Spec §4.2's prose and its figures disagree** about whether a tract can lose its last repeat.
3. **Spec §4.2 calls the part-repeat cutoff "10 base pairs"** where it counts a compressed rank.
4. **`arch/read_likelihoods.md` §4.2** sketches `stutter_rates_for` beside the distribution and
   names the context's field `truncated_mass_lost`; both moved.
5. **Spec §12's fourth test asks periods 1 to 6; the tripwire runs 2 to 6.**
6. **`Provenance`'s ladder puts *supplied* below *borrowed***, and `weaker_of` follows it.
