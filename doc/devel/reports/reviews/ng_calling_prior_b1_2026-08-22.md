# Code Review: ng_calling_prior_b1
**Date:** 2026-08-22
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B1 of the genotype-prior plan — the ported Dirichlet-multinomial primitive
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff of step B1 of
  [`calling_prior.md`](../../ng/impl_plan/calling_prior.md), branch `ng-calling-prior`. One file,
  +350/−5.
- **Reviewed against:** base commit `27eda128` plus `tmp/b1.patch`.
- **In-scope files:**
  [src/ng/calling/genotype_prior/dirichlet_multinomial.rs](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs)
- **Deliberately out of scope:** the `PriorRow` bundle it consumes (reviewed and fixed at A2);
  frozen production, read as reference and as an oracle.
- **Categories dispatched, and one of them is not a standard category.** `reliability`, `naming`
  and a combined `smells` + `idiomatic`. The fourth was written for this step: **port fidelity and
  numerics**, because the step's whole claim is that it computes what production computes, and no
  standing checklist asks that question. It returned the finding with the longest reach.

### 2. Verdict

**Request-changes.** One Blocker, four Major, sixteen Minor. Every Blocker- and Major-class finding
was demonstrated by a mutation or a measurement, not by reading.

### 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | `Finished dev profile … in 6.48s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `25 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `22 passed; 0 failed` |
| `cargo test --lib` | 0 | `4032 passed; 0 failed; 11 ignored` |

Findings labelled "Needs verification": **0.**

**Mutation totals across the three agents that ran them: 28 run, 7 survived, 4 changed no
behaviour** — reliability 7/4/1, port-fidelity 11/1/1, smells+idiomatic 10/2/2. Genuine survivors,
after each agent proved or disproved a behaviour change: **3**, all in the reliability report and
all rolled into the Blocker. The naming agent ran reproductions and compiled renames rather than
mutations.

### 4. Open questions and assumptions

1. **Which of production's two spellings should ng be bit-identical to?** They disagree by an ulp
   (B2 below). This review recommends staying with the shared primitive the plan names and
   *recording* the difference, because the alternative is a parity test against a private function
   in a frozen file. It is worth a ruling before the production differential is built. Affects
   **M1**.
2. **Should this module have a criterion bench?** The hot-path finding could not be settled either
   way because none covers `ng::calling::genotype_prior`. Affects the one Low-confidence finding,
   which proposes no code change.

### 5. Top 3 priorities

1. **B1** — every fixture runs at a reference concentration of 1, which is the corner the caller
   occupies only at one sample; three silently-wrong implementations pass the whole module.
2. **M1** — the provenance prose names the wrong production function and credits the port with a
   change production had already made.
3. **M2** — the port is bit-identical to the variant the STR path calls and differs by an ulp from
   the one the SNP/indel path runs, which is where the design's headline measurement came from.

### 6. Findings

#### Blocker

**B1: dirichlet_multinomial.rs:127 — every fixture in the file runs at a reference concentration of
1, and three wrong implementations pass because of it.** *Category: reliability.* Confidence: High.
All 25 tests draw their concentration from one helper, which pins the reference entry at 1 and the
table at six alleles or fewer. That is the leave-one-out concentration of **one sample at a
biallelic site**; what the primitive is actually handed is the run's seed plus the cohort's expected
allele copies (spec §6), which at a thousand diploid samples puts the reference entry near 2,000.
Three mutants pass the whole module: clamping every entry at 1, hard-coding the reference's entry,
and truncating the fold at six alleles. Measured, they move a row by **9.92 nats at a hundred
samples**, 16.71 at three thousand, and get **57 of 78 genotypes wrong** on a twelve-allele locus.
The recommended fix is to widen the **bit-parity** grid rather than the oracle's, since parity has
no tolerance to argue about; the agent verified parity holds at reference entries of 201, 801, 2001
and 6001, at ploidy 8 and 255, and out to twelve alleles.

#### Major

**M1: dirichlet_multinomial.rs:44–54 — the provenance names the wrong production function and
claims a change production had already made.** *Category: port fidelity.* The doc said the ported
primitive is "already shared between production's two callers" and that filling a caller slice is
what this port changed. Measured: `crate::genetics::dirichlet_multinomial_log_priors` has **one**
non-test caller, `src/ssr/cohort/em.rs:306`; the two hits in `posterior_engine.rs` are a doc comment
and a test. Production's SNP/indel engine runs its own copy, `fill_log_indep_per_g_from`, whose doc
says it computes the same thing "in place from `scratch.alpha` / `scratch.lgamma_alpha` to avoid a
per-record allocation" — **exactly the change this port claims as its own, for the same stated
reason.** The spec repeats the error at §9, so it originates upstream and is not this step's to fix.

**M2: dirichlet_multinomial.rs — the port matches one of production's two spellings and differs from
the other by an ulp, on the path the design's headline number came from.** *Category: port
fidelity.* `fill_log_indep_per_g_from` associates differently — it sums the per-allele terms and
adds the multinomial coefficient last, where the port folds from the coefficient. Replicating it and
comparing: **112 of 492 genotype values differ, by at most one unit in the last place** (largest
disagreement 7.1e-15 nats). So `the_port_matches_production_bit_for_bit` pins fidelity to the STR
cohort's variant, while the GIAB 83.6% → 94.6% measurement of spec §2.2 was taken on the other. It
cannot move a genotype; it will surprise whoever builds the production differential.

**M3: dirichlet_multinomial.rs — the oracle's `1e-12` tolerance is a fact about the grid, and the
doc says it is not slack.** *Categories: reliability, port fidelity (convergent).* Measured
disagreement between the two routes: 7.19e-14 at a reference entry of 201, 9.72e-13 at 801,
**2.05e-12 at 2001** and **7.62e-12 at 6001**. So the stated tolerance stops holding somewhere
between four hundred and a thousand diploid samples — harmless to genotyping, and ordinary `lgamma`
cancellation rather than a defect, but the claim as written is wrong.

**M4: dirichlet_multinomial.rs:249–253 — the doc says two tests guard the zero-count skip from
opposite sides, and only one does.** *Categories: naming, smells, port fidelity (three-way
convergent).* Run the mutation and `an_allele_a_genotype_does_not_carry_cannot_move_its_prior`
**passes**: its two concentrations perturb the affected rows by about `1e-22` against values near
`0.69`, six orders of magnitude below an ulp, so the two rows stay bit-identical with the branch
gone. The same test also passed a `+ α_a × 1e-18` perturbation planted on that exact branch. It
killed 2 of 8 behaviour-changing mutations and neither uniquely. **The skip is guarded by the parity
test alone** — which stops being able to guard it the day `src/genetics.rs` moves.

#### Minor

Sixteen, of which the ones that changed code: `concentration_of`'s alternatives sum to
`alt_total × allele_count / 2` rather than to `alt_total`, so at the grid's six-allele shape with
`0.5` they sum to 1.5, half again as much as the reference entry — the parameter is the *largest*
alternative, not a total, and it corrupts three assertion messages. The 2:1 test's doc promises to
fail "the moment anyone raises `α_ref`" while its own fixture hard-codes it three functions away,
and conflates two different ratios: at `α_ref = 10` the *marginalized* ratio is 19.998:1, and the
22:1 the spec also records is the plug-in path's Hardy–Weinberg ratio, a different quantity. The
monotonicity assertion in the hom-reference test cannot fail, because the `1 − 1.5θ ± 3θ²` windows
at the three θ do not overlap and the tolerance check already forces the ordering. "Tracks" is on
the project's banned-placeholder list. `checked`, `untouched`, `previous`, `oracle`, `row_for` and
`term` name what they hold poorly. A concentration entry of `0` or `+∞` is guarded only in debug and
becomes `−∞` or `NaN` in release with nothing recording it. And "transcription" reads as DNA→RNA to
this project's intended reader.

#### Nits

The row loop is three levels of nesting with two nested tuple patterns in fourteen lines, and the
order-of-operations promise the whole port rests on lives forty lines above the fold it constrains.
`pochhammer_ln` and `concentration_of` are unidiomatic names. The 1% band in the 2:1 test passes
with 1% of its tolerance to spare and breaks if a diversity is added to the loop.

### 7. Out of scope observations

- **Spec §9's reuse map has the same error as M1**, describing the primitive as shared between the
  engine and the STR path. The engine references it in prose and tests only.
- `itertools::izip!` would flatten the tuple nesting; itertools is not a dependency and adding one
  to remove two parentheses is not a trade worth making.
- **No criterion bench covers `ng::calling::genotype_prior`**, so the one hot-path question — that
  `lgamma(α_a + k_a)` is recomputed per genotype though only `alleles × ploidy` distinct values
  exist — has no evidence either way. It saves nothing at diploid biallelic (4 calls, 4 distinct
  values) and reaches about 5× at tetraploid with four alleles.

### 8. Missing tests to add now

1. **The parity grid at cohort-scale concentrations** — reference entries of 201, 2001 and 6001,
   ploidy to 8, alleles to 12 (B1).
2. **A skip test that does not read production** — hom-reference at a one-allele locus against
   hom-reference at a four-allele one, which the skip makes bit-identical and adding the zero terms
   back does not (M4).
3. **A reference entry off 1 in the 2:1 test**, so the ratio is shown moving when the pair moves.

### 9. What's good

- The arithmetic is a faithful transcription: **1,546,974 genotype values compared against
  production over ploidy 1–16 × 1–8 alleles × six concentrations, 0 bit mismatches.**
- The numerics are sound at the edges: the `MIN_ALT_CONCENTRATION` floor gives a finite row, ploidy
  255 and a 30-allele tetraploid table (40,920 genotypes) agree with the oracle to 7e-13 or better,
  and nothing overflowed or produced a non-finite entry anywhere in the committed range.
- The `3θ²` tolerance is genuinely tight — the true neglected term is `1.75θ²`, so the assertion
  uses 58% of its budget at each of the three diversities.
- Two oracles that catch different things: five mutations die only to bit-parity, two only to the
  rising-factorial oracle.
- The grid count assertion caught a number the author had guessed (1,264 against the real 348)
  before the review began.

### 10. Commands to re-verify

- `<worktree>/scripts/dev.sh cargo fmt --check`
- `<worktree>/scripts/dev.sh cargo clippy --lib --tests --all-features -- -D warnings`
- `<worktree>/scripts/dev.sh cargo test --lib`
- `<worktree>/scripts/dev.sh cargo test --release --lib ng::calling::genotype_prior`

Per-category files are left as an audit trail in `tmp/review_2026-08-22_ng-calling-prior-b1/`.
