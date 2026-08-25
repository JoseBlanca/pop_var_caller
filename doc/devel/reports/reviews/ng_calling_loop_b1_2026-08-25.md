# Code Review: ng calling loop — B1, the E-step for one sample

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step B1
**Implementation report:** [ng_calling_loop_b1_2026-08-25.md](../implementations/ng_calling_loop_b1_2026-08-25.md)
**Fixes applied:** [fixes_applied_2026-08-25_v3.md](fixes_applied_2026-08-25_v3.md)

## 1. Scope

The working-tree diff of B1 — three files, +644 lines: the new
`src/ng/calling/inference/summarise_condition.rs`, `SampleScoringBuffers` and its accessor in
`src/ng/calling/mod.rs`, and one `pub mod` line. `src/ng/calling/likelihood/` and
`src/ng/calling/allele_candidates/` are out of scope: two other branches own them.

**Five agents, each in its own git worktree, detached at `2b6aceae` with the patch applied.**
Categories: `reliability`; `errors` + `defaults`; `naming` + `idiomatic`; `smells` +
`refactor_safety` + `module_structure`; `extras` (hot path, diff-matches-intent) + the skill's
step 8a, the diff's own quantitative claims. **Merging nine checklists into five agents is a
deviation from the skill's one-agent-per-category rule**, taken because the diff is one small
new file; recorded rather than silent.

## 2. Verdict

**Request changes** — 3 Blockers, 7 Majors, 14 Minors, 11 Nits. Every Blocker and Major is
about what the tests *cannot see* or about a check that does not do what its doc says; none is
about the arithmetic, which every agent that checked it found correct.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --all -- --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | no warnings |
| `cargo test --lib` | 0 | `4534 passed; 0 failed; 14 ignored` (A2 left 4,528) |
| `cargo test --all-targets --all-features` | **red, pre-existing** | aborts at `benches/psp_writer_perf.rs:386`, the standing `PROJECT_STATUS.md` item |
| `cargo test --lib --bins --tests --examples --all-features` | **red, pre-existing** | `--example ng_generic_loci_dump` gives `2 passed; 11 failed`; **a detached worktree at `2b6aceae` without the patch prints the identical line** |
| `cargo test --release --lib ng::calling` | **red, pre-existing** | `493 passed; 4 failed` — four `should_panic` tests in `ng::calling::likelihood`, another branch's |

## 4. The three Blockers

**BL1 — `sample_scoring_buffers` picks two rows and nothing tests that it picks the right
ones.** Making the method ignore its `sample` argument entirely leaves every test in
`ng::calling` green. The test that looks like the one to catch it — the three-sample scoring
test, whose `assert_ne!` is commented *"sample 1 was scored on sample 0's likelihood row"* —
passes under the mutation, because it also redirects the expected-copies row, which
`score_one_sample` overwrites on every call, so each sample gets a different leave-one-out
term and therefore a different posterior **from the same reads**. This is the exact hazard the
row accessors' own doc says they exist to prevent.

**BL2 — the `posterior_row.len()` check is release-held, guards a real truncation, and no test
reaches it in either profile.** With it deleted, a three-genotype locus scored through a
two-entry posterior row returns copies `[1.333…, 0.667…]` where the answer is `[1.0, 1.0]` —
the reference allele over-counted by a third of a copy, no panic. **Both** invariants the
three-sample test asserts are satisfied by that wrong answer: the truncated row is renormalised
to one, and the copies sum to the ploidy because
`Σ_a Σ_g p_g·c_{g,a} = ploidy·Σ_g p_g` for *any* subset of genotypes whose posterior sums to
one.

**BL3 — `largest_score.is_finite()` cannot see a `NaN`, and the test that claims to reach that
case reaches a different one.** `largest_score` is only ever assigned through
`score > largest_score`, and every comparison against a `NaN` is false, so the maximum is never
itself a `NaN`: the check is in effect "the maximum is not `±∞`". The shipped test filled the
*whole* row with `NaN`, so what tripped it was the `−∞` the maximum started from — the panic
prints `came out -inf`. **A `NaN` in a non-maximal genotype passed.** Measured under
`--release`: likelihoods `[NaN, 0, 0]` returned normally with every posterior entry and every
expected copy `NaN`. In debug the only thing that stopped it was a `debug_assert!`.

**Three agents found BL3 independently** (`reliability`, `errors`, `extras`), which is why it
is the one finding here with no residual doubt.

## 5. The seven Majors

1. **A part-`NaN` likelihood row is not refused in release** — `errors`' statement of BL3, with
   the downstream trace: `SampleAlleleCopies::new`'s finiteness check is debug-only and
   `fill_sample_concentration`'s `max(0, ·)` returns the non-`NaN` operand, so the next pass's
   concentration comes back as the bare seed.
2. **`at_one_sample_the_concentration_comes_back_as_the_seed` cannot fail independently.**
   Replacing the leave-one-out term with `cohort − cohort` — which makes every sample's
   concentration the bare seed at **every** cohort size, so the loop is inert — leaves it
   green. It is killed only by the hand-computed test.
3. **No test runs a second pass**, so the one-sample fixed point that spec §13 test 1 names is
   not reached by anything.
4. **The scratch's `UNWRITTEN_SCRATCH_VALUE` sentinel does not survive to be seen here.** With
   a sample's own copies left unwritten, the release call returns normally and the
   concentration comes back as exactly the seed, because `f64::max` returns the other operand
   on a `NaN`.
5. **`pub` on `score_one_sample` is a dead-code silencer, not a visibility decision** —
   narrowing all three new items to `pub(crate)` produces exactly three "never used" errors,
   which is the only thing `pub` was buying.
6. **Step C1's flat prior-free pass cannot go through this function at all**, and in release it
   would silently become the seeded pass C1 exists to prevent. See §7.
7. **`SampleScoringBuffers` is the one bundle in this call chain with public fields and no
   constructor**, while `Concentration`, `CohortAlleleCopies`, `SampleAlleleCopies` and
   `PriorRow` all enforce their shape in `new` — and its own doc's claim that it is "made only
   by" the accessor is contradicted by a test in the same patch.

## 6. The diff's own numbers — 19 claims checked, 17 correct

The pattern is the inverse of the one this project's earlier milestones produced. **Every
number quoted about the author's own fixtures re-ran verbatim**, including all five mutation
strings, the hand-computed concentration, posterior and copies, the −691 floor
(`ln 1e-300 = −690.7755`), and the "does not compile" claim (two `E0499` plus one `E0502`).

**The two wrong ones are both prose about *properties*, not numbers about fixtures** — which is
worth recording, because the standing advice says to look at the author's own figures:

- *"the total is at least one … so the division cannot divide by zero"* — only if every score
  is finite, which the release-held check does not establish (BL3);
- *"A `NaN` or an `−∞` maximum"* — a `NaN` maximum is unreachable.

One claim is right and underspecified: *"about 430 reads at Phred 10"* is 434.3 for reads the
genotype does **not** explain, and 9,491 for reads it does.

## 7. The one finding that reaches past this step

**Step C1's flat pass.** Steps 1 and 2 of `score_one_sample` are unconditional, and spec §3
says the first pass of every round runs with **no prior at all**. Handing this function a flat
`GenotypePriorModel` does not achieve that: step 1 still reads `cohort_expected_copies`, which
on the first pass holds `UNWRITTEN_SCRATCH_VALUE`. Measured — it panics in debug
(`got [NaN, NaN]`) and in release `f64::max` swallows the `NaN`, so the concentration comes
back as the seed exactly, **which is the seed-only first pass spec §3 spends four paragraphs
ruling out**.

The `smells` agent implemented a repair — `PassPrior::{Flat, LeaveOneOut}` as a value rather
than a code path, one entry point and one `match`, `507 passed; 0 failed` — and it is a good
one. **It is not applied here**, because the plan puts the flat pass at C1 and the shape of
C1's entry point is C1's decision, not B1's. Carried as an open item, with the measurement
attached so C1 does not have to rediscover it.

## 8. Hot path — every figure measured, no finding

- **Zero allocations, proved from the linked binary** rather than asserted. The crate is
  `unsafe_code = "forbid"`, so a counting allocator cannot be written inside it; instead
  `score_one_sample`'s complete call inventory was read out of `objdump`: one
  `fill_sample_concentration`, one vtable `blr`, three `exp`, one `memset`
  (`sample_expected_copies.fill(0.0)`), and 15 cold panic paths. **No allocator symbol.** The
  shipped prior behind the seam is clean too.
- **No bounds checks in the hot loops** — no `panic_bounds_check` anywhere in the function;
  every loop is iterator-driven.
- **The `if copies != 0` guard earns its place**: 90 of the 126 copy-table entries are zero at
  six alleles, and removing it cost 8–29% of the E-step's own arithmetic across twelve runs
  against a stub prior. **Against the prior that actually ships it disappears** — 0–4%, inside
  the run-to-run drift, because the marginalized Dirichlet prior is about 250–280 ns of a
  380–405 ns call at six alleles. That is spec §2's *"the prior is not a small term … it is the
  part of the E-step that carries the expensive function"*, now with a number against it.
- **The release-held asserts cost nothing measurable** — the sign of the difference flipped
  between builds.

## 9. The range commitment

Nothing in the six tests went above 3 samples or 3 alleles. The `reliability` agent ran the
wide case rather than reasoning about it — **1,000 samples over 8 alleles, 36 genotypes, scores
spread over about 10,500 nats**, cohort copies the exact sum of the samples' own, one sample
carrying the rarest allele alone. Every posterior summed to 1 and every copies row to 2 within
`1e-12`. The count-path desync threshold is absolute at `−1e-6` while the rounding error in a
2,000-copy sum over 1,000 samples is about `1e-10`, four orders below it.

**So the ceiling hides nothing numeric**, and the agent's recommendation — spend the coverage
on the missing tests rather than on a bigger cohort — is what was done. What three samples does
hide is a *fixture-shape* problem: BL1 survives at three samples and would survive at a
thousand.

## 10. Out of scope

- `src/ng/calling/likelihood/` — four `should_panic` tests fail under `cargo test --release`,
  each targeting a `debug_assert`. Another branch's; it is why `--release` cannot be read as
  green module-wide, and it is the same blocker `PROJECT_STATUS.md` records against adding a
  release CI step.
- `src/ng/calling/genotype_prior/` — whether `SampleAlleleCopies::new` should refuse a
  non-finite entry in release is that module's owner's call. B1 works around it with a check of
  its own.
- `genotype_table.rs`'s `homozygous_alleles()` against `genotype_prior`'s
  `homozygous_allele_for` — one concept, two names, both predating this patch.

## 11. What's good

The arithmetic. Every agent that checked the four steps against spec §2 found them right, in
the right order, with the M-step correctly absent. The disassembly confirmed the design's two
hard requirements — no allocation, no bounds checks — rather than leaving them as intentions.
And the quantitative discipline held: 17 of 19 claims re-ran verbatim, with both failures in
prose about properties rather than in a number about a fixture.
