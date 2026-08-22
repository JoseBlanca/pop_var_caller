# Code Review: ng_calling_prior_b2
**Date:** 2026-08-22
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B2 of the genotype-prior plan — the two-branch inbreeding mixture
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `d163c6d8` on branch `ng-calling-prior`, +430/−2 across two files.
  Step B2 of [`calling_prior.md`](../../ng/impl_plan/calling_prior.md). **The commit was made
  without a review**, at the owner's direction to close the previous session; this review is the
  one it owed.
- **Reviewed against:** `d163c6d8`, with `c12a5d45` (step B1) as the base.
- **In-scope files:**
  [dirichlet_multinomial.rs](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs)
  (the mixture, the log-sum-exp helper, seven tests) and
  [mod.rs](../../../../src/ng/calling/genotype_prior/mod.rs) (`PriorRow::ploidy`, the re-export).
- **Deliberately out of scope:** step B1's primitive and its four tests (reviewed and fixed at B1);
  the rest of `mod.rs` (reviewed at A2), except where B2 made one of its premises load-bearing;
  `src/genetics.rs`, `src/var_calling/`, `src/ssr/` — frozen production, read as reference and as
  oracles.
- **Categories dispatched:** `reliability` (always); `naming` (always, and this project's rules are
  stricter than the checklist); `smells` + `idiomatic` (the no-allocation contract); `errors` +
  `defaults` + `refactor_safety` (the module has no `Result` by design, and this is the seam's
  first implementation); and a fifth written for this step — **model fidelity and numerics**,
  because B2's headline is a deliberate departure from what it ports and no standing checklist asks
  whether a mathematical claim is true.

### 2. Verdict

**Request-changes.** Two Blockers (one issue, found independently by two agents), six Major, six
Minor. **The step's mathematics survived the attack**; everything filed is about what the code
leaves untested and what its prose claims wrongly.

### 3. Execution status

Run by the orchestrator, in the container, on `d163c6d8`:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 4.33s` |
| `cargo test --lib` | 0 | `4039 passed; 0 failed; 11 ignored; finished in 710.16s` |

Reproduced by the agents in their own worktrees: `32 passed` debug and `29 passed` release for
`ng::calling::genotype_prior`, matching the commit.

Findings labelled "Needs verification": **0.**

**Mutation totals across the five agents: 44 run, 22 survived, 15 changed no behaviour.** Genuine
behaviour-changing survivors after each agent proved or disproved the change: **3**, all folded
into B1 and M2. Per agent — reliability 11/6/4, model fidelity 8/3/2, smells+idiomatic 13/9/4,
errors+defaults 6/2/1; the naming agent compiled every rename it proposed instead of mutating.

### 4. Open questions and assumptions

1. **Should a genotype the prior rules out carry `−∞` or the probability floor?** The code writes
   `−∞`; spec §8, spec §12 test 3, arch §1.1 and this plan's own step B2 line all say floor.
   Production's mixture writes `−∞` too, so the code matches what it ports and the four documents
   describe `wright_genotype_log_priors` rather than the engine. It is not cosmetic: the comparator
   at step F1 is ported from that Wright function, which floors, so the two implementations behind
   one seam would differ by convention as well as by model. **Owner's ruling.** Affects **M3**.
2. **Is `PriorRow::ploidy` right to return a bare `u32`** where `GenotypeTable` and
   `GenotypeTableView` both return the checked `Ploidy` newtype? Two agents split on it — the
   newtype compiles and passes but adds two `expect` paths to an accessor on a per-sample
   per-pass path. Affects **M5**, whose fix is taken from the cheaper side.

### 5. Top 3 priorities

1. **B1** — the seam's only implementation has no test that goes through the seam; an
   implementation that dropped the inbreeding coefficient entirely left the whole module green.
2. **M1** — every mixture fixture but two pins the reference concentration at 1.0, the single value
   at which a wrong `Σα` is invisible. This is B1's own Blocker, repeated in the same file.
3. **M2** — a mutation the implementation report records as killed is not killed, and the doc
   comment that rests on it gives a mechanism that cannot occur.

### 6. Findings

#### Blocker

**B1: dirichlet_multinomial.rs:186 — the seam's only implementation is reached by no test**
**Categories:** reliability, errors+defaults, smells+idiomatic (cross-category). **Confidence:** High.

`MarginalizedDirichletPrior`, its `GenotypePriorModel` impl and the `mod.rs` re-export are the only
items B2 makes reachable from outside the file, and all seven new tests drive the private
`fill_marginalized_log_priors` instead. `mod.rs`'s two trait tests use a local stand-in.

Two agents independently replaced the impl body with a call to the random-mating primitive —
dropping the coefficient — and the module stayed green at 38 passed. The mutant is not
behaviour-preserving: at `F = 0.95`, `α = (1, 0.01)`, biallelic diploid, the real row gives a
het:hom-alt ratio of **0.051** and the mutant **1.980**, a factor of 39, about 16 on the Phred
scale. At `F = 1` through the trait the heterozygote goes from impossible to ordinary.

This is where every future caller enters, and where step F1's second implementation will be
compared against this one, so a defect here is invisible to the whole comparison rather than to one
test.

**Fix:** a test that builds the row through `MarginalizedDirichletPrior` and `InbreedingF`, pinned
bit-for-bit against the bare-coefficient spelling, plus one assertion that the coefficient reached
the mixture at all.

#### Major

**M1: dirichlet_multinomial.rs:641-668 — every mixture fixture but two pins the reference
concentration at 1.0.** **Categories:** reliability, model fidelity. **Confidence:** High.

The concentration this function is handed is the leave-one-out one (spec §6): near 1 at one sample,
near 2,000 at a thousand diploid ones. The normalisation identity — the test the commit message
calls "what the two branches being on one scale *is*" — runs nine shapes × four `F` × three
diversities, all at a reference entry of 1.0. Measured with a `Σα` that hard-codes that entry: the
identity's worst error is **4.0e-15** against its `1e-9` budget, so it cannot see the defect at all;
raise the reference to 6,001 and the same identity misses by **0.95**. The mutant was killed only by
the two Wright tests, which use a large total for an unrelated reason.

**This is the defect B1's review already found once**, in the same file — B1's parity grid was
widened from a pinned 1.0 to `{1, 201, 2001, 6001}` for exactly this reason, and B2's new fixtures
went back to the narrow value.

**Fix:** add the reference entry to the normalisation grid.

**M2: dirichlet_multinomial.rs:238-254 — a guard documented as load-bearing that no test can
defend, and a mutation wrongly recorded as killed.** **Categories:** reliability, model fidelity,
naming, smells+idiomatic (all four). **Confidence:** High.

The doc says "**Both short-circuits are load-bearing rather than defensive** … without the first,
`−∞ − −∞` would make every row `NaN`", and the implementation report §5 records "the `−∞`
short-circuit dropped from `log_sum_exp_2` | killed — every row becomes `NaN` at `F = 0`". Measured
independently by four agents: delete the first, the second, or both, and every row entry is
bit-identical and the module stays green (40 passed debug / 37 release with both gone). `−∞ − −∞`
needs *both* arguments infinite, which takes `F = 0` and `F = 1` at once or a concentration entry
of zero that `Concentration::new` refuses.

The sentence is also backwards about which guard covers which end: at `F = 0` it is the **second**
that fires, not the first.

**Fix:** correct the doc to what the guards do buy — the both-`−∞` pair returning `−∞` instead of
`NaN`, and a saving of two `exp` and a `ln` per homozygote at `F = 0` — add a direct unit test of
the helper, and correct the implementation report's mutation table.

**M3: dirichlet_multinomial.rs:217 — `−∞` reaches the caller, against four design documents.**
**Categories:** errors+defaults, reliability (cross-category). **Confidence:** High.
**See open question 1 — this is the one finding whose resolution is the owner's.**

At `F = 1` the heterozygote's entry is `f64::NEG_INFINITY`, and it is reachable through the
**public** seam today, not only through the bare-coefficient door: `InbreedingF::try_new(1.0)`
returns `Ok`. Spec §8, spec §12 test 3, arch §1.1 and this plan's step B2 line all call for
`PROBABILITY_FLOOR`. Production's own mixture writes `−∞` as well, so the code matches what it
ports; what the four documents describe is `wright_genotype_log_priors`, whose floor is baked in and
which is the comparator's source at F1.

**Fix applied for now:** state in the trait's contract what the implementation does and that the
question is open, so a second implementer meets it there rather than in one implementation's tests.
The ruling itself is deferred to the owner.

**M4: dirichlet_multinomial.rs:199 — the mixture checks nothing about the coefficient, in either
profile.** **Categories:** errors+defaults. **Confidence:** High.

`fill_marginalized_log_priors` takes a bare `f64` and goes straight to `inbreeding.ln()` and
`(1.0 - inbreeding).ln()`. Fed real values, identically in debug and release: `NaN` gives a wholly
`NaN` row; `1.5` the same; and **`-0.1` gives a row that is `NaN` on both homozygotes and finite on
the heterozygote** — a half-poisoned row, which normalises to a plausible wrong answer. Nothing
panics.

The newtype holds the line today, so the defect is latent. It stops being latent at F1, where a
second implementation in this same file will reach for the same bare-float helper. And the doc's
justification — that this is a test-only path — is not true: the trait implementation routes every
caller through it.

**Fix:** one `debug_assert!` that the coefficient is a fraction in `[0, 1]`, which is where A2 puts
every check on a *value*; the inclusive upper bound stays, because admitting `1` is the function's
stated purpose.

**M5: mod.rs:312 — `ploidy()`'s premise is checked nowhere, and B2 is what made it load-bearing.**
**Categories:** errors+defaults. **Confidence:** High.

The doc says the first genotype's counts "cannot disagree with the table it came from". They can:
`PriorRow::new` takes `genotype_allele_counts` as a bare slice with no tie to any table and checks
only its length. Measured on a diploid biallelic table whose first row was edited to `[6, 0]`:
`ploidy()` returned 6 and the homozygous-reference entry moved **5.89 nats**, with nothing raised in
either profile. With the first row zeroed, `ploidy()` returned 0 — a value `Ploidy::try_new` refuses
— and the scale correction silently vanished, reproducing exactly the row that the dropped-correction
mutation produces and that three tests catch when it arrives as *code*.

B2 made this matter: the correction sits on the identical-by-descent branch only, so a wrong `m` is
not a shared constant the row may carry — it re-weights the mixture.

**Fix:** a release `assert!` that the first genotype's counts do not sum to zero (that one is
silent, so it earns release), and a debug `debug_assert!` that every genotype agrees on the total.

**M6: the implementation report's account of production is wrong in its citation and understates
the size by three orders of magnitude.** **Categories:** model fidelity. **Confidence:** High.

The report says `pipeline.rs:343` "passes the cohort's *fitted* coefficient". That line passes the
**CLI knob**, whose default is 0; the fitted per-sample values arrive at the following lines through
`with_fixation_index_overrides`. The conclusion holds and is stronger than written.

The size is understated because production feeds the mixture the leave-one-out concentration, so
`Σα` grows with the cohort and the inflation factor `Σα(Σα + 1)` grows with its square. Measured on
the same fixture the report uses — biallelic diploid, tomato1's fitted diversity of 6 in 10,000,
het:hom-alt ratio at `F = 0.8`:

| | outbred | correct at F = 0.8 | uncorrected |
|---|---|---|---|
| 1 sample | 2.00 | 0.222 | 0.400 |
| 50 samples | 188.7 | 0.493 | 181.8 |
| 1,000 samples | 1818 | 0.499 | 1816.5 |

At one sample the coefficient still does 90% of its work; at 50 samples **3.6%**, at 1,000 samples
**0.09%**. So in production the fitted coefficient is not weakened, it is very nearly inert, and the
report's "about 1.8 times too likely" is the *smallest* case rather than a representative one.

#### Minor

- **Mi1: dirichlet_multinomial.rs:727-748** — `the_wright_agreement_closes_as_the_concentration_grows`
  seeds its accumulator with `f64::INFINITY`, so the first of its four totals asserts
  `gap < INFINITY`. Three of four comparisons were real. *(reliability)*
- **Mi2: mod.rs:312-316** — `PriorRow::ploidy()` is a new public accessor with no direct test; a
  mutant reading the last genotype instead of the first survives (correctly — it is a no-op — but
  nothing pins which row is read or the slice bound). *(reliability)*
- **Mi3: dirichlet_multinomial.rs:184** — `MarginalizedDirichletPrior` derives nothing, where every
  other public type in the module derives `Copy, Clone, PartialEq, Debug`, so a
  `&dyn GenotypePriorModel` cannot be printed or recorded. The seam exists for a two-way comparison;
  a result that cannot name the prior that produced it is not auditable. *(errors+defaults)*
- **Mi4: naming** — `fill_marginalized_log_priors` names what it shares with its sibling rather than
  what differs: both fill functions marginalize, and what this one adds is the inbreeding mixture,
  which is what its own doc, its helper and five of seven tests call it. Not pinned in any design
  document. *(naming)*
- **Mi5: naming** — `log_outbreeding` / `log_inbreeding` for `ln(1 − F)` and `ln F` introduce a
  second vocabulary for the pair that three lines below is called `independent_draws` /
  `identical_by_descent`, which is the spec's §3.2 wording. "Outbreeding" also usually names
  heterozygote excess — a negative `F` — in the literature. *(naming)*
- **Mi6: naming** — one quantity under three names: `scale_of_the_random_branch` in the code,
  `shared_constant` in a test, "a shared additive constant" in the prose. The code's name says where
  the value goes, not that it is the row's offset from a true log-probability. *(naming)*

#### Nits

The `1e-9` mass tolerance has about 250,000× headroom on its own fixtures and the test's doc gives
no measured figure, where B1's tolerances all do. `mod.rs:76` splits the `pub mod` block with a
`pub use`, which the plan asks be kept clean for a sibling branch's merge. The six-buffer
`PriorRow::new` setup is written out three times in the test module. `log_sum_exp_2`'s `a.max(b)`
shift is unreachable at the concentrations this caller commits to — the worst `|b − a|` is 86.9 nats
against the 709 `exp` needs — and a `larger = a` mutant survives with last-place differences only.
The `ln_1p` two-term spelling of the helper would use one `exp` instead of two, but would break the
port's fidelity to production's spelling. `concentration[usize::from(allele.0)]` panics with a bare
"index out of bounds" where every other refusal in the module names both buffers it compared.

### 7. Out of scope observations

- **`InbreedingF::try_new` accepts `1.0`**, which spec §7 says it should not. The tightening to
  `[0, 1)` belongs to [`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md)
  Milestone A. B2's doc claims the tightening as though it had happened, which is corrected here.
- **Production's mixture has the scale defect and it is live.** Outside the frozen tree's edit
  permission and outside this plan, but it should not be left only in a commit message. Now
  measured at cohort scale under **M6**.
- **The spec is unamended.** §3.1's "the constant cancels" and §3.2's mixture still describe the
  uncorrected form; the code doc says so and the spec edit is the owner's.

### 8. Missing tests to add now

All eight were added and are green; each was re-run against the mutant it targets.

| test | catches |
|---|---|
| `the_seam_and_the_bare_coefficient_agree_and_both_carry_the_inbreeding_coefficient` | B1 |
| `the_seam_rules_out_heterozygotes_at_the_greatest_coefficient_the_newtype_accepts` | B1, and signals when `InbreedingF` is tightened |
| `log_sum_exp_2_returns_the_finite_argument_when_the_other_is_impossible` | M2 |
| `a_haploid_row_is_unmoved_by_the_inbreeding_coefficient` | a wrong ploidy or `Σα` in the correction, at the shape where both branches must coincide entry by entry |
| `a_monomorphic_locus_is_unmoved_by_the_inbreeding_coefficient` | the same, at one allele |
| `ploidy_returns_the_copy_count_every_genotype_sums_to` | Mi2 |
| `a_first_genotype_carrying_no_copies_is_refused` | M5, in release |
| `genotypes_that_disagree_on_the_copy_count_are_refused_in_debug` | M5, in debug |

### 9. What's good

- **The Wright oracle earned its place.** It is the only check in the plan that exercises both
  branches at once, and it is what found the scale defect that four documents and production all
  carry. Everything else in B1 and B2 passes either way.
- **The lying-lookup test** — handing a real diploid table a lookup that says nothing is homozygous
  — is the right way to pin "one function decides homozygosity", and it fails the obvious inline
  spelling rather than merely describing it.
- **`at_no_inbreeding_the_mixture_leaves_the_random_mating_row_untouched` asserts bit equality**, not
  a tolerance, which is what makes "the correction cannot move an outbred sample" checkable.
- **The correction is added to the identical-by-descent branch rather than subtracted from the
  other**, which is why the `F = 0` row is bit-identical to the primitive's and nothing B1 pins
  moved.
- **The implementation report volunteered the missing review** and named the finding as something to
  attack rather than confirm. Four of the five agents attacked it; the mathematics held.

### 10. Commands to re-verify

```
scripts/dev.sh cargo fmt --check
scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
scripts/dev.sh cargo test --lib ng::calling::genotype_prior
scripts/dev.sh cargo test --release --lib ng::calling::genotype_prior
scripts/dev.sh cargo test --lib
```

The release run is not optional: three of the module's checks are held in release deliberately and a
debug run cannot tell `assert!` from `debug_assert!`.
