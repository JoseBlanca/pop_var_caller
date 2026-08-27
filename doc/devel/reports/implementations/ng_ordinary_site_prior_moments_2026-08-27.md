# ng — the ordinary-site prior's two numbers: implementation report

**Branch** `ng-prior-moments`, worktree `../pop_var_caller-prior-moments`, cut from `main` at
`9f15f5e5`.

**Plan** [`../../ng/impl_plan/ordinary_site_prior_moments.md`](../../ng/impl_plan/ordinary_site_prior_moments.md).
**Design authority** [`../../ng/spec/ordinary_site_prior_moments.md`](../../ng/spec/ordinary_site_prior_moments.md),
with [`../../ng/spec/ordinary_site_seed.md`](../../ng/spec/ordinary_site_seed.md) §3 for the
identity that turns two moments into a concentration pair.
**The measurements this work rests on** were made before it started and are in
[`../ng_ordinary_site_prior_moments_2026-08-27.md`](../ng_ordinary_site_prior_moments_2026-08-27.md);
no step here is a sweep.

**One report, one section a step.** The plan's steps are small and several are deletions, so a
file apiece would be a file of two paragraphs; each section below carries the step's own contract,
what shipped, and what was measured about it.

---

## A1 — the population's mean alternative-allele frequency, in closed form

**Contract (plan A1).** Add `p_fixed_alt + p_segregating · a/(a+b)` beside the heterozygosity
integral that already exists, with the same population-not-panel framing. Its own commit, because
a wrong mean frequency is a plausible number at every panel size and nothing downstream refuses it.

### What shipped

`FrequencyDensity::expected_alternative_frequency`, in
[`src/ng/parameter_estimation/joint/fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs),
one line of arithmetic beside `expected_heterozygosity`. Two of the density's three parts
contribute: positions where the population carries only a non-reference base are at frequency one
and carry their whole share, positions that segregate contribute their Beta's mean, and positions
carrying only the reference base contribute nothing.

Nothing calls it yet. A2 is where the seed starts reading it.

### The oracle: the search's own answer, where the search is exact

**The cheapest available proof that the replacement is right is the thing it replaces.** At one
diploid individual the panel has three allele-count classes, two of them free once normalised,
against the two-parameter family's two parameters — so the search reproduces those classes exactly
and its mean frequency is the density's own (spec §9, report §9's third reading). At larger panels
it is fitted over more classes than it has parameters and drifts, up to 1.22× the truth at 200
individuals, so the comparison is available only here.

Measured, on six densities at one individual and no inbreeding, the search's mean frequency over
the closed form's:

| density | search | closed form | ratio |
|---|---:|---:|---:|
| tomato-like, `Beta(0.20, 1.00)` | 1.664544e-3 | 1.666667e-3 | 0.9987 |
| human-like, `Beta(0.35, 1.20)` | 1.459505e-3 | 1.461290e-3 | 0.9988 |
| flat, `Beta(1.00, 1.00)` | 2.999262e-3 | 3.000000e-3 | 0.9998 |
| the unit tests' lopsided fixture, `Beta(0.50, 2.00)` | 2.798057e-2 | 2.800000e-2 | 0.9993 |
| middling, `Beta(4.00, 4.00)` | 2.999262e-3 | 3.000000e-3 | 0.9998 |
| reference base rare, `Beta(3.00, 0.60)` | 4.337823e-3 | 4.333333e-3 | 1.0010 |

**The band is 0.9987× to 1.0010×**, an order of magnitude inside the search's own 1% resolution
(`SearchPrecision::fast`). The two routes share no algebra: one maximises a log-likelihood over
Beta-binomial class weights, the other is one line of Beta moments.

**This test dies at A5**, with the search it is measured against. What survives as the permanent
check is the hand-computed one below.

### What the fixtures share, and what was added because of it

**The five densities this repository already sweeps all have `a ≤ b`** — the alternative allele
rare, or the two shapes balanced. On every one of them, reading `b/(a+b)` for `a/(a+b)` returns a
number too *high*, so a reader checking the sign of an error would see a consistent story and
conclude the formula was merely mis-scaled. `Beta(3, 0.6)` — the population where the reference
base is the rare one at the positions that vary (report §2) — was added for that reason, and the
swap moves the answer **down** there and up on the other five.

**Six rows, five distinct answers.** `Beta(1, 1)` and `Beta(4, 4)` are both symmetric, so both have
mean a half and both give 3.000 in 1,000. They differ in spread and not in mean, so for this
quantity they are one fixture. The set is kept whole because it is the set the earlier measurements
used; the duplication is recorded rather than left to be found.

### Mutations run

Four, on the formula itself, each against all three of the step's tests
(`the_expected_alternative_frequency_is_the_densitys_own`,
`the_two_point_masses_carry_their_own_ends`,
`the_closed_form_frequency_is_the_searchs_own_answer_at_one_individual`):

| mutation | tests failing |
|---|---|
| `a` and `b` swapped | 2 of 3 |
| the fixed-non-reference share dropped | 3 of 3 |
| the invariant share read in its place | 3 of 3 |
| the segregating share dropped | 3 of 3 |

The swap leaves `the_two_point_masses_carry_their_own_ends` green, which is correct: at
`p_fixed_alt = 1` and at `p_invariant = 1` the Beta contributes nothing either way, so that test
pins the two ends and says nothing about the shape between them.

The source was restored from a backup and the restore checked with `git diff` before anything
else ran.

### Validation

All in the container, from this worktree:

- `cargo fmt --all -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,877 passed, 0 failed, 14 ignored**, against the branch's 4,874 at
  `9f15f5e5`. Three tests added, none removed.
- `cargo doc --no-deps --lib` — **27 unresolved links, the same 27 as at `9f15f5e5`**. The crate
  denies broken intra-doc links, so this command exits 101 on the pre-existing set; what matters is
  that the count did not move.
