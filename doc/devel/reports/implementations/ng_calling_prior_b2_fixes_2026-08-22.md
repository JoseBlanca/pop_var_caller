# ng genotype prior — B2 review fixes: eight tests, three corrected claims, one question for the owner

*Fix-application report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Applies
[the B2 review](../reviews/ng_calling_prior_b2_2026-08-22.md) to step B2 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md).*

## 1. What the review found, in one line

**The mathematics held; the tests and the prose did not.** The corrected mixture is the exact
mixture up to one genotype-independent constant at every ploidy and allele count — checked against
a rising-factorial oracle over 1,080 rows, worst spread 2.5e-11 nats — and the `F = 0` row is
bit-identical to the primitive's over 100,220 entries, far beyond the six shapes the step tested.
Everything applied below is about what the step left unpinned, and about three of its own claims
that were measurably wrong.

## 2. Findings and their disposition

| id | finding | disposition |
|---|---|---|
| B1 | the seam's only implementation is reached by no test | **Applied** |
| M1 | every mixture fixture but two pins the reference concentration at 1.0 | **Applied** |
| M2 | a guard documented as load-bearing that no test can defend; a mutation wrongly recorded as killed | **Applied** |
| M3 | `−∞` reaches the caller, against four design documents | **Applied** — owner ruled at Checkpoint B: floor it |
| M4 | the mixture checks nothing about the coefficient | **Applied** |
| M5 | `ploidy()`'s premise is checked nowhere | **Applied** |
| M6 | the production account: wrong line, size understated | **Applied** (report only) |
| Mi1 | the closing test's first iteration cannot fail | **Applied** |
| Mi2 | `PriorRow::ploidy()` has no direct test | **Applied** |
| Mi3 | `MarginalizedDirichletPrior` derives nothing | **Applied** |
| Mi4 | `fill_marginalized_log_priors` names what it shares, not what differs | **Applied** |
| Mi5 | `log_outbreeding` / `log_inbreeding` are a second vocabulary for one pair of branches | **Applied** |
| Mi6 | one quantity under three names | **Applied** |
| — | `PriorRow::ploidy()` should return the `Ploidy` newtype | **Deferred** — see §5 |
| — | a `name()` method on the trait so a run can record which prior produced a row | **Deferred to F1** — see §5 |
| — | the test module's three copies of the six-buffer setup | **Deferred** — see §5 |
| — | `InbreedingF` tightened to `[0, 1)` | **Out of scope** — prerequisites plan, Milestone A |

## 3. What changed

**Eight tests added.** The module goes from 32 passed debug / 29 release to **40 / 36**.

| test | the defect it catches |
|---|---|
| `the_seam_and_the_bare_coefficient_agree_and_both_carry_the_inbreeding_coefficient` | B1 — bit equality between the trait path and the bare-coefficient path at four shapes and four coefficients, plus one assertion that the coefficient is used and not merely passed |
| `the_seam_rules_out_heterozygotes_at_the_greatest_coefficient_the_newtype_accepts` | B1 at the model's edge. It also acts as a tripwire: when the prerequisites plan tightens `InbreedingF`, this constructor stops returning `Ok` and the test's `expect` says so |
| `log_sum_exp_2_returns_the_finite_argument_when_the_other_is_impossible` | M2 — the helper's guards, which no row fixture can fail on |
| `a_haploid_row_is_unmoved_by_the_inbreeding_coefficient` | one copy is one draw, so both branches are the same statement and `F` must be inert. **Every haploid genotype is homozygous**, so this is the only test in which the identical-by-descent branch fires on every entry — the two branches have to coincide entry by entry, not merely add to one |
| `a_monomorphic_locus_is_unmoved_by_the_inbreeding_coefficient` | the same at one allele, where `α_0 / Σα` is exactly 1 |
| `ploidy_returns_the_copy_count_every_genotype_sums_to` | Mi2 — the slice bound at every shape including the one-allele table, and the premise it rests on |
| `a_first_genotype_carrying_no_copies_is_refused` | M5 in release |
| `genotypes_that_disagree_on_the_copy_count_are_refused_in_debug` | M5 in debug |

**Two checks added to `PriorRow::new`** (M5). A first genotype whose copies sum to zero is refused
**in release**, because that one is silent: `ploidy()` returns 0, the correction becomes
`lgamma(Σα) − lgamma(Σα)` — exactly zero — and the mixture reverts to the unscaled one this step
exists to fix. Genotypes that disagree on the total are refused **in debug**, which is where this
module puts every check on a *value*.

**One check added to the mixture** (M4): a `debug_assert!` that the coefficient is a fraction in
`[0, 1]`. Measured before it existed, identically in debug and release: `NaN` and `1.5` each give a
wholly `NaN` row, and `-0.1` gives a row that is `NaN` on both homozygotes and finite on the
heterozygote — a half-poisoned row, which normalises to a plausible wrong answer.

**The normalisation grid now runs the reference entry to 6,001** (M1), which is a thousand diploid
samples' worth of leave-one-out counts. Worst error over the whole grid is 1.5e-11 against the
unchanged `1e-9` budget.

**The closing test is seeded with a measured bound** rather than `f64::INFINITY` (Mi1), so all four
of its totals carry an assertion. The seed is `1e-1`: the correct code's gap at `Σα = 100` is
1.10e-2 and the uncorrected mixture's is **7.98e-1**, so the seed is what makes a dropped correction
fail on the first iteration rather than the second — verified by running that mutation, which now
fails with *"at Σα 100 the gap to Wright was 0.798…"*.

**Three renames** (Mi4–Mi6): `fill_marginalized_log_priors` → `fill_inbreeding_mixture_log_priors`;
`log_outbreeding` / `log_inbreeding` → `log_weight_independent_draws` /
`log_weight_identical_by_descent`, matching the branch names three lines below;
`scale_of_the_random_branch` → `shared_normalising_term`, which is what the test and the prose
already call it.

**Four doc claims corrected**, each because it was measured false, not because it read badly:

- `log_sum_exp_2`'s "both short-circuits are load-bearing … `−∞ − −∞` would make every row `NaN`".
  Deleting either guard or both leaves every entry bit-identical. The doc now says what they do buy
  and notes which guard fires at which end — at `F = 0` it is the **second**, which the original had
  backwards.
- "matches the Wright formulas to one part in a million". Measured **1.9e-5** over the test's own
  grid, which is why the tolerance is `1e-4`.
- `ploidy()`'s "cannot disagree with the table it came from". It can, and now it is checked.
- The mixture's size claim, replaced with a measured table across cohort sizes — see §4.

**One new measured paragraph** on where `lgamma(Σα + m) − lgamma(Σα)` stops being computable: it is
a difference of two nearly equal numbers of order `Σα·ln Σα`. Measured as the row's departure from
one unit of probability, biallelic, ploidy 2 and 8, `F` to 0.95: 1.5e-11 at `Σα = 7.2e3` — past the
top of the committed cohort range — 9.1e-11 at 1.2e5, 1.1e-9 at 1.2e6, 2.8e-7 at 1.2e8. Nothing in
range is affected; the 1.2e6 figure is why the normalisation identity is not also run at the total
the Wright oracle uses.

## 4. The production defect is three orders of magnitude larger than the step reported

The step reported it as making a heterozygote "about 1.8 times" too likely. **That is the one-sample
figure, and it is the mildest case.** The defect inflates the random-mating branch by
`Σα(Σα + 1)`, and the concentration production feeds the mixture is the leave-one-out one, which
grows with the cohort. Measured on the step's own fixture — biallelic diploid, tomato1's fitted
diversity of 6 in 10,000, heterozygote-to-homozygous-alternative prior ratio at `F = 0.8`:

| | outbred (`F = 0`) | correct at `F = 0.8` | uncorrected | how far the coefficient got |
|---|---|---|---|---|
| 1 sample | 2.00 | 0.222 | 0.400 | 90% |
| 50 samples | 188.7 | 0.493 | 181.8 | 3.6% |
| 1,000 samples | 1818 | 0.499 | 1816.5 | 0.09% |

So on a cohort **production's fitted inbreeding coefficient is very nearly inert.** The report and
`PROJECT_STATUS` are corrected; production itself is frozen and out of this plan.

Also corrected: the report cited `pipeline.rs:343` as the line passing the fitted coefficient. That
line passes the CLI knob, whose default is 0. The fitted per-sample values arrive at the lines
after it, through `with_fixation_index_overrides`.

## 5. Deferred, and why

- **`PriorRow::ploidy()` returning the `Ploidy` newtype.** Two review agents split on it: the
  newtype compiles and passes, but it adds two `expect` paths to an `#[inline]` accessor on a
  per-sample per-pass path. The premise is now checked at construction, which is where the cost is
  paid once rather than per call, so the newtype buys type-safety at the seam's boundary and nothing
  else. Worth revisiting when the calling loop exists and the real call frequency is known.
- **A `name()` method on `GenotypePriorModel`.** The seam exists for a two-way comparison and a
  result should name the prior that produced it. But adding a required method to the trait now
  ripples into F1's implementation and into the loop, which is beyond this step. `Debug` is derived
  in the meantime, so a `&dyn GenotypePriorModel` can at least be printed. **Recorded for F1.**
- **The test module's three copies of the six-buffer setup.** A real duplication, and the reviewer
  compiled an extraction that works. Minimal-diff discipline: this step's fixes are already the
  largest change to the file since it was written, and a test-helper refactor landing in the same
  commit would make the behaviour-relevant parts harder to read. **Recorded for whoever next
  touches the module's tests.**

## 6. M3, and the ruling — the floor goes in

**A genotype the prior rules out was written as `−∞`, and four design documents asked for the
probability floor.** Spec §8, spec §12 test 3, arch §1.1 and this plan's own step B2 line all say
floor. **Production's mixture writes `−∞`** (`safe_ln` returns `NEG_INFINITY`), so ng matched what
it ports — what the four documents describe is `wright_genotype_log_priors`, whose floor is baked in
and whose only shipping caller is the paralog filter.

**Owner's ruling, 2026-08-22: floor it.** The reason it is not cosmetic is step F1:
`PlugInWrightPrior` is ported from that Wright function and will floor. Two implementations behind
one seam would otherwise differ by convention as well as by model, and the seam exists to attribute
a genotype difference to the model.

**The floor is asymmetric, on purpose.** `1 − F` is floored; `F` is not.

- `1 − F` can reach a **row entry**: at `F = 1` it is a heterozygote's only branch, so an unfloored
  `ln 0` writes `−∞` into the output. Floored, the entry is about −694.7 and every entry of the row
  is finite.
- `F`'s own `ln 0` at `F = 0` never reaches a row entry — it makes the identical-by-descent branch
  impossible, and `log_sum_exp_2`'s second short-circuit returns the other branch exactly. Flooring
  it too would send **every outbred sample** down the general path, paying two `exp` and a `ln` per
  homozygous genotype to move nothing. Measured: flooring both leaves the outbred row bit-identical,
  so the cost would buy literally zero change.

**The outbred row is untouched by the floor**, which is what makes this safe: zero bits of
difference across 72 combinations of ploidy, allele count, diversity and reference concentration.
The floor sits 300 orders of magnitude below anything the row carries. At `F = 1` the two
homozygotes keep their concentration ratio to 4 parts in 10¹⁶.

**Correction to this report's first draft**, which recommended *against* flooring on the grounds
that it would break the bit-exact outbred test. Measured: it does not. That claim was reasoned
rather than run, and it was wrong.

**Also corrected: `F` can be exactly 1 today.** The two `F = 1` tests are not testing an
unreachable edge — `InbreedingF::try_new(1.0)` returns `Ok` on this commit, so the shipping seam
admits it. Once the prerequisites plan tightens the newtype to `[0, 1)`, the largest coefficient it
can carry is `1 − 2⁻⁵³`, where `1 − F` is about `1.1e-16` and the floor never bites — at which point
the floor becomes a guard for the bare-coefficient function, which admits `1` by design so the limit
stays testable.

The trait's contract now states the rule for **every** implementation rather than describing one.

## 7. Mutations re-run after the fixes

Every one applied against the fixed tree, in the container, and every one now dies:

| mutation | before | after |
|---|---|---|
| the seam ignores the inbreeding coefficient | survived (38 passed) | **killed** — both seam tests |
| `Σα` hard-codes the reference entry at 1.0 | killed only by the two Wright tests | **killed** — and now by the normalisation identity too |
| both `−∞` short-circuits deleted | survived (40 passed) | **killed** — the direct helper test |
| the scale correction dropped | killed by 3 tests | **killed** by 4, and the closing test now fails on its first total rather than its second |
| the ploidy hard-coded at 2 | killed by 1 test | **killed** by 3, including at haploid and monomorphic shapes |
| the identical-by-descent branch always reads allele 0 | killed by 3 | **killed** by 4 |
| the floor on `1 − F` removed | — (added by the M3 ruling) | **killed** — both `F = 1` tests |

The two new `PriorRow::new` checks were confirmed to fire in the profile each claims: a zero-summing
first row is refused in debug **and** release; disagreeing rows are refused in debug and not in
release.

## 8. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.82s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `40 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `36 passed; 0 failed` |
| `cargo test --lib` | 0 | see the commit message |

The pre-existing red gates named in the plan's session brief (`--all-targets` clippy, the
`psp_writer_perf` bench panic, 17 unresolved intra-doc links) are untouched and live in files this
branch does not open.
