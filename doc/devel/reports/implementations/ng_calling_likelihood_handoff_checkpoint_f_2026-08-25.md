# ng read likelihoods — handoff at Checkpoint F

*2026-08-25. Branch `ng-calling-likelihoods`, worktree `../pop_var_caller-calling-likelihoods`,
plan [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). Written for
whoever picks this up next. It follows
[the Checkpoint C/D handoff](ng_calling_likelihood_handoff_checkpoint_cd_2026-08-24.md), which is
now a historical record — **its eight things still hold**, and this adds what Milestones E and F
cost.*

## Where the plan stands

**Twelve of fifteen steps are done: A1–A2, B1–B2, C1–C2, D1, E1–E3, F1–F3.** That is the whole
generic path, the whole stutter distribution, and the whole STR emission seam with two models in
it. **Three remain: G1, H1, H2** — the censored term and the STR row.

| milestone | steps | state |
|---|---|---|
| A–D — the generic path | A1…D1 | ✅ (Checkpoint C/D) |
| E — the stutter distribution's three changes | E1, E2, E3 | ✅ |
| F — the STR emission seam | F1, F2, F3 | ✅ **Checkpoint F reached** |
| G — the censored term | G1 | ☐ next |
| H — the STR row | H1, H2 | ☐ |

**Everything through F3 is merged to `main`** (fast-forward, `081f81d2`); `105cba0c` is F3 on top.

Suite: **4,392 passed, 0 failed, 14 ignored** in the library target; **182 in
`ng::calling::likelihood`**. Build and test with
`/Users/jose/devel/pop_var_caller-calling-likelihoods/scripts/dev.sh cargo test`; validate with
`cargo clippy --lib --all-features --tests -- -D warnings` — **not** `--all-targets`, which is red
on `main` in `examples/ng_duplicated_class_harness.rs` and `benches/freebayes_bookkeeping.rs`.

## The nine things a successor must not rediscover

**1. On the STR path the locus *is* the tract.** An observation's bases and a candidate's bases
are the repeat run **alone, without flanks** — the generator slices exactly that
(`region_seq[tract]`, [`locus_generation/ssr.rs`](../../../../src/ng/locus_generation/ssr.rs)), and
the owner confirmed it on 2026-08-25. **This is the opposite of the generic path**, where an allele
is the whole locus as a carrier has it, and F1 shipped with the generic reading in its doc before
F2 caught it. Getting it wrong compares flanks against flanks and puts the tract at the wrong
offset for every resize.

**2. Two places encode the same contraction boundary, and they must move together.** A tract must
keep **at least one repeat**: `StutterModel::unreachable_mass` counts a total contraction as
unplaceable, and `enumerate_placements` refuses a placement that would empty the tract. **Change
one without the other and the report silently stops describing the scoring.** The boundary is not
arbitrary — it is the only reading of spec §4.2 that reproduces the two sizes §4.2 itself states
(2.0 parts in a million and 2.0 parts in a thousand); the other reading gives 1.0 in ten million
and 1.0 in a thousand. Milestone E's review re-derived both independently and confirmed it.

**3. `unreachable_mass` understates the loss on interrupted tracts, and that is unfixed.** It takes
a period and a **total** repeat count; reachability is **per run**. A tract of two runs of two
repeats cannot lose three even though it holds four, and the report counts that contraction as
reachable. Closing it needs the run structure to reach the distribution. **This is the sharpest
open question in the plan** — it is a number the row uses to keep candidates comparable.

**4. Model A and Model B disagree on part-repeat lengths by about ninefold, and that is correct.**
A charges the fitted part-repeat share; B has no such branch and absorbs the odd base as end slop
at the flat per-base rate. On a ten-base read A prefers the four-repeat candidate (1.869e-4 against
1.094e-4) and B the three-repeat one (9.706e-4 against 1.167e-5). **Do not "fix" it** — the
part-repeat branch is why spec §4.1's comparison chose Model A. Where the length difference *is* a
whole number of repeats the two agree completely, and that is the check.

**5. Two names moved and the architecture still uses the old ones.**
`truncated_mass_lost` → **`unreachable_mass`** (only the cutoffs are truncation, and at period 1
nothing is truncated at all — the part-repeat branch does not exist there, and that term is 2 in
100 against 1 in 10¹³ for a cutoff tail). `ssr_emission.rs` → **`stutter_rates.rs`** for the
parameter adapter; the real `ssr_emission.rs` is F1/F2/F3's. `arch/read_likelihoods.md` §4.2 still
says `truncated_mass_lost` and still sketches `stutter_rates_for` beside the distribution, which is
the placement E3 rejected on module-edge grounds.

**6. `alignment` and `parameter_estimation` are siblings that import each other nowhere**, checked
in both directions across the whole crate including test-only code. That is why the fit-to-rates
adapter lives under `calling::likelihood` and not beside the distribution: `StutterModel`'s own
contract says it fits nothing.

**7. The review fan-out earns its cost, again.** Over E and F it ran 26 mutations and found: a
zero-repeat candidate answered as though it held one (the largest mis-scaling that function can
produce, arriving silently — now a `NonZeroU32`, so a compile error); a test asserting only
`>= 0.0` where every wrong implementation passes; a doc naming a guard that could not fail on its
behalf; and the fit's own ceiling untested, where the model floors and the loss reports exactly
zero for a row that has over-allocated six parts in a hundred. **None of those crash.**

**8. Check the review's own suggested fixes.** One proposed `assert_eq!(lost, 0.0)` "because the
clamp fires at every cell". It fires at **12 of 18** — at period 1 with a single repeat the model
can still place only 0.51, so 0.49 is a real loss rather than a negative one to clamp away.
Verifying the suggestion is what found the more interesting fact.

**9. Two operational things that cost time.** The full suite has hit a **transient rustc internal
compiler error** in the codegen backend twice, on unrelated examples (`ng_catalog_window_probe`,
`ng_inbreeding_resolution`); it did not reproduce either time — **re-run before diagnosing**. And
**always drive the build by the worktree's absolute path**: the shell's working directory drifted
to the main checkout once mid-session, and an edit landed there instead. Recovered, but check
`git -C /Users/jose/devel/pop_var_caller status` if anything looks odd.

## What G1 needs, so it is not read twice

Spec §5.2, and the plan's own bullet: the **factorised** form on pure candidates — the letter match
on the witnessed prefix times the closed-form tail `P(length ≥ ℓ | a)`, both geometric tails capped
at E2's cutoffs — and the **exact sum** over reachable stretchings on interrupted candidates.

Three tests, and the first is the one to get right:

- **the complement identity** — `P(≥ ℓ) + P(< ℓ)` equals the truncated distribution's own total,
  which is **one minus `unreachable_mass`, not 1** (spec §12 test 12). Item 2 and item 3 above are
  both live here.
- where the constraint admits exactly one length change, censored equals complete **bit for bit**;
- **a partial never out-discriminates a complete observation** on a stated parameter set (test 13).

`SsrEmissionModel::censored_emission` is `unimplemented!()` today, in both Model A and the oracle.
That is deliberate: a stub returning a number would let a truncated read out-discriminate a whole
one, which is the exact failure test 13 exists to catch.

## The six questions waiting for the owner

None blocks G or H.

1. **`unreachable_mass` understates the loss on interrupted tracts** (item 3).
2. **Spec §4.2's prose and its figures disagree** about whether a tract can lose its last repeat.
   The branch follows the figures, in two places now (item 2).
3. **Spec §4.2 calls the part-repeat cutoff "10 base pairs"** where it counts a compressed rank;
   ten of those admit about 13 base pairs at period 4.
4. **`arch/read_likelihoods.md` §4.2** carries two stale names and one rejected placement (item 5).
5. **Spec §12's fourth test asks periods 1 to 6; the sums-to-one tripwire runs 2 to 6**, because at
   period 1 the part-repeat branch is unreachable by construction and the total is *supposed* to
   fall short.
6. **`Provenance`'s ladder puts *supplied* below *borrowed***, and the new `weaker_of` follows it.

## The standing habit, and this run's evidence for it

**Every number about your own work is measured before it is written.** Wrong-first-time on this
run: a sweep's widest loss guessed at 1.1e-3 that measured 2.06e-2; a short-tract loss guessed at
2.5e-3 that was 2.525e-3; a clamp said to fire at 14 cells that fires at 12; a placement count of 1
that is 2; a test-count baseline of 4,358 that was 4,364; and a doc claiming the sums-to-one
tripwire catches a dropped complement, which it cannot. **Every one was caught by asserting it in a
test rather than stating it in prose.** Do that.
