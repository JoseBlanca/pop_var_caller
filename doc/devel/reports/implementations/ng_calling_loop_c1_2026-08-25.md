# ng calling loop — C1: the flat first pass

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step C1
**Design authority:** [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md) §3, §7, §13
test 3; [arch/calling_em_loop.md](../../ng/arch/calling_em_loop.md) §2.1, §4
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Two category reviews raised eight Majors.
> **One was a defect**: the flat arm bypassed two release-held shape checks, so a mis-shaped
> buffer that made the seeded arm panic made the flat arm return a wrong posterior in silence.
> **Three were claims in this report that were wrong** — the window's cohort-size axis, the
> rejected design's failure mode, and, largest, §4's account of the trap: **it is a delay of
> eight passes, not a different answer**, and `spec/calling_em_loop.md` §3's stronger claim was
> not reproduced on any fixture tried.
> [The review](../reviews/ng_calling_loop_c1_2026-08-25.md), and
> [what was done about it](../reviews/fixes_applied_2026-08-25_v5.md).

---

## 1. Plan

**The first pass through a locus runs on the reads alone.** The reason is mechanical: the
leave-one-out prior is built from the cohort's expected allele copies, and those are what a
*previous* pass produces. On the first pass there is none.

**The choice of *how* to say that was settled by B1's review**, which measured what the obvious
alternative does. Handing the seam a flat `GenotypePriorModel` does **not** give a prior-free
pass: step 1 still runs and reads a buffer holding the `NaN` sentinel — a panic in debug, and in
release a concentration equal to the bare seed, because `f64::max` returns the other operand on
a `NaN`. That is precisely the seeded first pass §3 exists to prevent. So the choice is a
**value**, matching `arch/calling_em_loop.md` §2.1's own rule that every switch of the design is
a value rather than a code path.

## 2. Changes made

- **`PassPrior<'a>`** — `Flat`, or `LeaveOneOut { model, inbreeding }`. `score_one_sample` takes
  it in place of a `&dyn GenotypePriorModel` and an `InbreedingF`.
- **Steps 1 and 2 became a `match`.** The flat arm writes a zero prior row and builds no
  concentration, so nothing reads the cohort's expected copies. It writes the row rather than
  branching inside step 3, so the hot loop keeps one spelling; the cost is `genotypes` stores
  against `genotypes` calls to `exp`.
- **The own-copies finiteness check is now the seeded arm's**, gated by one predicate,
  `PassPrior::reads_the_cohort`. On a flat pass that row is *expected* to hold the sentinel —
  that is the whole situation the flat pass exists for — and step 4 overwrites it without ever
  reading it.

## 3. Deviations

**One, and it is in the test helper rather than the code.** `run_passes` took the starting
cohort as the *sum* of the samples' own copies. That made the counterfactual untestable: with a
cohort of `n` samples' worth, every sample's leave-one-out term is `n − 1` samples' worth on
pass 1, which is not what a seeded first pass is. Spec §3 defines the alternative as *"the seed
concentration on its own — the prior with its cohort term set to zero"*, so the helper now takes
the starting cohort explicitly, and the trap test sets it equal to one sample's own copies.

**This was found by the fixture failing**, not by reading: the first version of the trap test
had the seeded arm calling heterozygous too.

## 4. Tests added

Three.

| test | what it pins |
|---|---|
| `the_flat_first_pass_needs_no_cohort_summary` | the flat pass runs at a locus whose cohort summary is still the `NaN` sentinel — the situation it exists for |
| `a_seeded_first_pass_converges_to_no_variant_where_a_flat_one_finds_it` | spec §13 test 3's **trap**: `0/0` against `0/1` at every sample, at 3 samples and at 63 |
| `after_the_flat_pass_the_copies_reflect_the_reads_and_not_the_seed` | test 3's first half — the copies the second pass's prior is built from |

**The counterfactual mutation is the oracle.** Replacing the flat arm with the seeded one — build
the concentration, run `MarginalizedDirichletPrior` — fails **exactly these three tests and
nothing else**: `23 passed; 3 failed`.

**The review added four more**, and renamed the trap test to what it actually shows: the flat
pass touching neither the cohort summary nor the concentration (the only release-profile cover
for the property the variant exists for), the two shape checks the flat arm had been bypassing,
and a silent sample's flat-pass contribution.

### The numbers behind the fixture, measured rather than assumed

- **The prior gap under the seed `[1.0, 0.000_5]` is 7.6009 nats.** Read off the row:
  `[0.6931471805599453, -6.907755278982137, -7.600402584500431]`.
- **The tests use a likelihood advantage of 1 nat**, about 4.3 Phred — one alternative read
  among a handful.
- **The effect has a window with two axes** *(corrected after review — this report first gave
  only one)*. Swept over advantage × cohort size (3, 6, 20, 63): at **2 nats and above the two
  starts agree at every size**. At **0.5 nats** they part at 20 and 63 samples but **not** at 3
  or 6, where the flat start does not reach the heterozygote either. At **1 nat** they part at
  every size tried.
- **And the difference is a delay, not a different answer** *(the review's largest finding)*. At
  63 samples the seeded start flips to heterozygous at **pass 9**, and both starts reach
  0.767332 copies per sample by pass 30; at three samples the flip is between passes 10 and 16.
  A rare-variant fixture — carriers among 60 firmly homozygous-reference samples — showed the two
  agreeing in every cell. **`spec/calling_em_loop.md` §3's claim that the seeded loop "converges
  to no-variant" was not reproduced**, and is raised for the owner.

**What survives unaltered is the mechanical reason for the flat pass** — the cohort's copies do
not exist on the first pass — which is what makes `PassPrior` necessary whatever §3's stronger
claim turns out to be worth. The delay is not nothing: production's comment records
expectation-maximization converging in 3 to 5 passes, so a locus needing 9 under one start and 1
under the other is a locus whose answer depends on where the cap falls. Spec §12's question 7 is
the measurement on real data.

## 5. Validation

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4603 passed; 0 failed; 14 ignored`. Before C1: **4,596**, so C1 adds
  seven: three as written and four the review added.
- `cargo test --release --lib ng::calling --all-features` — `557 passed; 0 failed`.

## 6. Trade-offs and follow-ups

- **The flat arm leaves `sample_concentration` and `prior_per_allele_workspace` untouched**, so
  they hold whatever the previous locus left. The review confirmed neither can reach an answer —
  `prepare_for_locus` refills both per locus, `fill_sample_concentration` rewrites the
  concentration whole, and both shipped priors write the workspace before reading it — and added
  the test that asserts the concentration is untouched.
- **⚑ A sample with no reads votes for a 50% allele frequency on the flat pass**, and neither §3
  nor §7 says whether it should. §7's *"the prior decides it alone"* has no prior to appeal to
  here. Measured at 63 samples: three silent samples put 3.0027 alternative copies into pass 1
  against 0.0030 with none. At three reads a position roughly one sample in twenty is silent at
  any position. Pinned by a test, raised for the owner.
- **At one sample the flat pass has nothing to do** — spec §7 — because the cohort term is zero
  at every pass, so a flat pass 1 followed by a seeded pass 2 reaches what a single seeded pass
  would have. It costs one pass and changes no genotype. Whether to spend a branch on skipping
  it is spec §12's question 6, not settled here.
- **`PassPrior::Flat` has no caller outside tests.** D1's loop is the caller, and it is what
  makes the "runs at the start of every outer round" half of §3 real; C1 builds the value, not
  the schedule.
