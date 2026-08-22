# ng genotype prior — B2: the inbreeding mixture, and a scale the spec licensed dropping

*Implementation report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step B2 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone B.*

## 1. Plan

`MarginalizedDirichletPrior`: spec §3.2's two-branch inbreeding mixture over B1's primitive, and
the first implementation of the step-8 seam. With probability `F` a sample's two copies are one
ancestral copy counted twice, so the genotype is homozygous for allele `a` with probability
`α_a / Σα`; with probability `1 − F` they are independent draws and the genotype is worth what the
Dirichlet-multinomial gave it.

Design authority: [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §3.2 (the mixture),
§3.3 (what homozygous means above diploidy, and the one-function rule), §7 and §12 test 3 (the
`F = 1` limit), §12 tests 1 and 2 (carried through the seam);
[`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §3.2.

## 2. The finding, and the one deliberate departure from what was ported

**The two branches were not on the same scale, and the spec is what licensed it.**

§3.1 has the primitive drop `lgamma(Σα + m) − lgamma(Σα)` because it is the same for every
genotype and cancels when the loop rescales the row. That is true of a row on its own. It is
**false the moment a second branch is mixed in**: the identical-by-descent term `α_a / Σα` is a
true probability, so mixing it with a row carrying that offset inflates the random-mating branch
by `Σα(Σα + 1)` at diploidy — and the inbreeding coefficient does a fraction of the work it
should.

**Measured three ways before anything was changed.**

- Subtracting the term from the primitive's row makes it sum to one, so the offset is exactly what
  separates the two scales. *(The four digit-strings originally quoted here could not be reproduced
  by the review and are withdrawn rather than restated. The identity is now asserted by a test
  instead of stated in prose —
  `the_mixed_row_is_a_probability_distribution_once_its_shared_constant_is_removed`, whose worst
  measured error is 1.5e-11 against a 1e-9 budget, over concentration totals from about 1 to
  7,200.)*
- Against the Wright formulas, the as-ported mixture is wrong by whole nats across
  `p ∈ {0.05, 0.2, 0.5}` and `F ∈ {0.25, 0.5, 0.8}`, and the error **grows with `F`** while the
  value itself barely moves with it — the signature of the inbreeding branch being swamped.
  *(The "0.29 to 2.20 nats" range originally given here mixes two different measures — the
  heterozygote contrast alone and the normalised row — and understates the worst single genotype.
  Withdrawn; the claim it supported is carried by the Wright oracle, which is a test.)*
  Corrected, the same comparison lands within **1.9e-5** at `Σα = 1e6` and **1.1e-2** at
  `Σα = 1e2` — the finite-concentration gap that closes as the limit is approached, and the reason
  the oracle's tolerance is `1e-4` rather than the "one part in a million" this report first
  claimed.
- **At the concentration this caller ships** — one sample, tomato1's fitted diversity of 6 in
  10,000, biallelic diploid — the het-to-hom-alt prior ratio comes out at `0.400` where the model
  says `0.222` at `F = 0.8`, and `0.200` against `0.105` at `F = 0.9`. A heterozygote made about
  **1.8 times** as likely as it should be, roughly 2.6 on the Phred scale, in the direction this
  caller is already weakest.

**Production has the same defect, and it is live rather than latent.** `posterior_engine.rs` mixes
`log_indep_per_g`, which carries the offset, with `log_p_effective[a] = ln(α_a) − ln(Σα)`, which
does not. It has gone unnoticed because `DEFAULT_INBREEDING_COEFFICIENT` is `0`, where the branch
short-circuits away and nothing shows — but the pipeline also hands the engine the per-sample
coefficients the diversity estimator *fitted*, through `with_fixation_index_overrides`, and those
are not zero on an inbred panel.

> **Corrected by the B2 review, 2026-08-22.** This paragraph first cited `pipeline.rs:343` as the
> line passing the fitted coefficient. That line passes the **CLI knob**, whose default is 0; the
> fitted per-sample values arrive at the lines immediately after it, as overrides.

> **The size was understated, and by three orders of magnitude.** See §9 below.

**Corrected in ng, owner-authorised 2026-08-22.** The correction is two `lgamma` calls per row,
added to the identical-by-descent branch rather than subtracted from the other — the same mixture
up to the shared constant a row is allowed to carry, and it leaves an outbred sample's row
bit-identical to the primitive's, so nothing B1 pins moves. **B1's own bit-parity with production
is untouched**: the primitive did not change, only what the mixture does with it.

**What found it was the oracle the plan asked for.** The Wright formulas at `F = 0.5` are the only
check in the plan that exercises both branches at once. Everything else in B1 and B2 passes either
way, because at `F = 0` the two spellings are identical.

**Owed to the spec:** §3.1's "the constant cancels" needs the qualification that it cancels in a
row and not in a mixture, and §3.2 should say the random-mating branch enters on the probability
scale. A spec edit is the owner's.

## 3. Changes made

| file | + | − |
|---|---|---|
| `src/ng/calling/genotype_prior/dirichlet_multinomial.rs` | 428 | 2 |
| `src/ng/calling/genotype_prior/mod.rs` | 15 | 0 |

- **`MarginalizedDirichletPrior`**, the first implementation of `GenotypePriorModel`, re-exported
  from `genotype_prior`.
- **`fill_marginalized_log_priors(row, inbreeding: f64)`** — the mixture on a bare coefficient.
  It exists so `F = 1` can be tested after the prerequisites plan tightens `InbreedingF` to
  `[0, 1)`, at which point no newtype will be able to carry a 1. *(Renamed to
  `fill_inbreeding_mixture_log_priors` by the review: both fill functions marginalize, and what
  this one adds is the mixture. The tightening has not happened —
  `InbreedingF::try_new(1.0)` returns `Ok` today — so `F = 1` is reachable through the trait as
  well, and both spellings are now pinned.)*
- **`log_sum_exp_2`**, ported from the same engine, with both `−∞` short-circuits. ~~They are the
  ordinary cases rather than edge ones: at `F = 0` the identical-by-descent branch is `−∞` on
  every homozygous genotype, at `F = 1` the other branch is `−∞` everywhere, and without the first
  short-circuit `−∞ − −∞` would make every row `NaN`.~~ *(The last clause is wrong — see §5. At
  `F = 0` it is the second short-circuit that fires, and deleting either or both leaves every row
  bit-identical. What they buy is the both-`−∞` pair returning `−∞` instead of `NaN`, and a saving
  of two `exp` and a `ln` per homozygote.)*
- **`PriorRow::ploidy()`** — read off the first genotype's copy counts rather than stored, since
  every genotype's counts sum to the ploidy. The mixture needs it for the correction's `m`.
  *(This report claimed the first row "cannot disagree with the table it came from". It can: the
  counts arrive as a bare slice with no tie to any table. A first row edited to `[6, 0]` moved a
  log-prior 5.89 nats with nothing raised. `PriorRow::new` now refuses a zero-summing first row in
  release and disagreeing rows in debug.)*

## 4. Tests added

Seven.

| test | what it pins |
|---|---|
| `the_mixed_row_is_a_probability_distribution_once_its_shared_constant_is_removed` | **The identity the correction *is*.** The two branches' weights add to one, so the row minus its shared constant sums to one — at nine shapes from haploid monomorphic to octoploid, four inbreeding coefficients and three diversities. The only non-diploid test here, and it had to be added: a correction hard-coding the ploidy at 2 passed everything else. |
| `the_concentrated_limit_matches_the_wright_formulas` | **The second independent oracle** (spec §3.2). At a concentration scaled to dominance the Dirichlet-multinomial collapses to a binomial draw, and the mixture must then say what Wright's biallelic diploid formulas say. Run at `F = 0` and `F = 0.5`, three frequencies, compared as differences between genotypes since the row carries a shared constant. |
| `the_wright_agreement_closes_as_the_concentration_grows` | That the `1e-4` above is the limit's own error rather than a tolerance chosen to pass: each tenfold rise in `Σα` must shrink the gap at least fivefold, at four totals a decade apart. *(As written, three of the four comparisons were real — the accumulator started at `f64::INFINITY`, so the first asserted `gap < INFINITY`. It is now seeded with a measured bound and all four carry an assertion.)* |
| `at_no_inbreeding_the_mixture_leaves_the_random_mating_row_untouched` | Bit-for-bit at six shapes and four diversities — so everything B1 pins survives the seam, and the correction is provably not shifting an outbred sample. |
| `at_full_inbreeding_only_the_homozygotes_survive_and_they_stand_at_the_concentration_ratio` | **Spec §12 test 3.** The heterozygote reaches `−∞` and the two homozygotes stand at `α_ref : α_alt`. |
| `the_homozygous_branch_reads_the_tables_lookup_rather_than_the_copy_counts` | **Spec §3.3's one-function rule, checked by lying to it.** A real diploid biallelic table is handed a lookup saying nothing is homozygous; both homozygotes must then come back with the independent-draws branch alone. An implementation comparing copy counts against the ploidy — the obvious inline spelling — would ignore the lie and fail. |
| `raising_the_inbreeding_coefficient_moves_mass_onto_the_homozygotes` | Monotone across five coefficients, and the two ends named: about 2:1 for an outbred sample, under 0.2 at `F = 0.95`. |

## 5. Mutation results

**Four run, four killed.**

| mutation | outcome |
|---|---|
| the scale correction removed (production's behaviour) | killed — the Wright oracle, its convergence check, and the normalisation identity |
| the ploidy hard-coded at 2 in the correction | **survived until the normalisation identity was added**; killed now |
| the homozygous branch reading the copy counts instead of the lookup | killed — the lying-lookup test |
| the `−∞` short-circuit dropped from `log_sum_exp_2` | ~~killed — every row becomes `NaN` at `F = 0`~~ **wrong; see below** |

> **Corrected by the B2 review, 2026-08-22 — this row was not killed, and the reason given for it
> cannot occur.** Four review agents independently deleted the first short-circuit, the second, and
> both: every entry of every fixture stayed bit-identical and the module stayed green. `−∞ − −∞`
> needs *both* arguments infinite, which takes `F = 0` and `F = 1` at once, or a concentration
> entry of zero that `Concentration::new` refuses. The sentence was also backwards about which
> guard covers which end — at `F = 0` it is the **second** that fires. So this step's mutation
> score was **three killed and one survivor**, not four killed, and `log_sum_exp_2` now has a
> direct unit test because no row fixture can fail on it.

## 6. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | `Finished dev profile … in 6.57s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `test result: ok. 32 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `test result: ok. 29 passed; 0 failed` |
| `cargo test --lib` | 0 | `test result: ok. 4039 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 627.50s` |

## 7. Review status

**Reviewed 2026-08-22**, by a five-agent fan-out that was given the scale-correction finding above
as something to attack rather than confirm. **The mathematics held** — the corrected mixture is the
exact mixture up to one genotype-independent constant at every ploidy and allele count, checked
against a rising-factorial oracle over 1,080 rows, worst spread 2.5e-11 nats. What the review found
was what the step left untested and what this report claimed wrongly:
[the review](../reviews/ng_calling_prior_b2_2026-08-22.md),
[what the fixes changed](ng_calling_prior_b2_fixes_2026-08-22.md).

## 8. Trade-offs and follow-ups

- **ng's mixture now differs from production's deliberately**, which is the first such divergence
  in this plan. A production differential on an inbred panel will show it, and should.
- **The spec owes two sentences** (§2 above).
- **`PROBABILITY_FLOOR` is still unimported.** The `F = 1` heterozygote reaches `−∞` rather than a
  floor, which after rescaling is a weight of exactly zero — an impossible genotype rather than a
  very unlikely one, which is the wanted answer. Its consumer is the comparator at F1.
- **The plug-in comparator will need the same correction** if it is to be compared fairly, since
  it too mixes a Hardy–Weinberg row with the identical-by-descent branch. Noted for F1.

## 9. The size of the production defect, measured at cohort scale

*Added by the B2 review, 2026-08-22. §2 above gave this as "about 1.8 times" a heterozygote's true
likelihood. That is the **one-sample** figure and the smallest one; the defect grows with the
square of the concentration total, and the concentration production feeds the mixture is the
leave-one-out one, which grows with the cohort.*

Same fixture as §2 — biallelic diploid, tomato1's fitted diversity of 6 in 10,000 — reading the
heterozygote-to-homozygous-alternative prior ratio, which is what the coefficient is for. An
outbred sample sits at 2:1 and an inbred one should fall well below it.

| | outbred (`F = 0`) | correct at `F = 0.8` | uncorrected at `F = 0.8` | how far the coefficient got |
|---|---|---|---|---|
| 1 sample | 2.00 | 0.222 | 0.400 | 90% |
| 50 samples | 188.7 | 0.493 | 181.8 | 3.6% |
| 1,000 samples | 1818 | 0.499 | 1816.5 | 0.09% |

So on a cohort **production's fitted inbreeding coefficient is very nearly inert**: at fifty samples
it moves the ratio 3.6% of the way to where the model puts it, and at a thousand, 0.09%. The 1.8×
of §2 is the corner where the defect is mildest.
