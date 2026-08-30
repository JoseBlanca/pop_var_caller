# ng calling loop — F1: the SNP/indel parity oracle

**Step:** F1 of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — ng's frequency loop
against the shipping caller's.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §10 (the
oracle and what it is allowed to excuse) and §13's eighth test.
**Date:** 2026-08-27. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

**ng's loop is a port, and a port has to prove it agrees with what it was ported from.** Given the
same genotype log-likelihood table, the same prior concentration, the same inbreeding coefficient
and the same allele count, the shipping caller's expectation-maximization loop and ng's now call
the same genotypes — asserted over nine fixtures in `src/ng/calling/loop_parity.rs`. **Two
differences exist and neither is excused**: each was found by this oracle failing, each is traced
to a decision already written down, and each is now asserted in the direction its record predicts,
so losing either correction breaks a test rather than quietly restoring agreement.

## 2. What is held identical, and why none of it is a transcription

**The likelihood table is the fixture's**, handed to both sides. Neither computes it, so a
difference between the two emissions cannot be mistaken for a difference between the two loops —
which is what §10 asks for in so many words.

**The concentrations are the same construction rather than two numbers that happen to agree.**
Production derives `α = [ALPHA_REF, θ/k, …]` from its `nucleotide_diversity`; ng derives it from a
`SpectrumSeed`'s reference concentration and alternative total, floored the same way. Seeding ng at
`(ALPHA_REF, θ)` makes the two arrays equal entry for entry, and **both read the same `ALPHA_REF`
out of `crate::genetics`** — so a change to production's constant moves both sides and the test
keeps meaning what it says. A fixture that typed the pair in would have been comparing two
transcriptions.

**Both sides are reached through code their own pipelines run.** Production through
`run_em_columnar`, which is what `run_em_for_record` calls; ng through `run_frequency_loop` and
`summarise_final_pass`, which is what its driver calls. The driver's own `call_locus` is *not* used,
and that is the point: it would build the likelihood table from the evidence, and the table is the
input under test.

**It needed no edit to the frozen tree.** `run_em_columnar`, `EmInputs`, `MergedAllelesView`,
`RecordScratch` and `RecordScratch::empty` are already `pub(crate)`; `PosteriorEngineConfig`,
`RecordLocus`, `EmDiagnostics`, `AlleleSupportStats`, `MergedAllele` and `mod backends` are `pub`.
The one earlier parity module cost production a `mod shape` declaration; this one costs nothing, and
the dependency still runs one way — ng reads production, production names nothing in ng.

**One thing is not independent, and it is worth knowing.** Production's winning genotype *index* is
resolved to allele copies through **ng's** genotype table, because production's own `GenotypeShape`
is not reachable from here. That is safe only because `genotype_table_parity` pins the two
enumerations value for value — but it does mean this file alone could not detect an
enumeration-order slip, and the sibling module is what does.

## 3. The two differences, both found by failing

**The inbreeding mixture — the one place the port departs on purpose.** Production mixes its
random-mating and identical-by-descent branches on two different scales: the first is a log-prior up
to a shared additive constant, the second a true probability, so the random-mating branch is
inflated by `Σα(Σα + 1)` and the inbreeding coefficient does a fraction of the work it should. ng
adds the missing `lgamma(Σα + m) − lgamma(Σα)` and it does all of it (owner, 2026-08-22; the
measurement is in `MarginalizedDirichletPrior`'s own documentation).

**At `F = 0` the branch short-circuits away on both sides**, which is why every outbred fixture
agrees exactly and why production's own default of zero hides this. At `F = 0.9` over ten samples
the oracle found it: a sample whose reads put its heterozygote **1.0 nat** above its homozygous
reference is called `0/0` by ng and `0/1` by production. The test asserts both calls by name and
also asserts that the coefficient moves ng at all — without that last check the fixture could show
a departure that was really a prior doing nothing on either side.

**The pass count — a difference in what is counted, not in where either stopped.** Both loops begin
with one E-step on the reads alone, before any prior: production's `EmStepPhase::FirstIteration`,
ng's initialisation. Production counts it as iteration 1; ng does not, because its `passes` counts
the passes that had a prior. **Measured on the fixture that takes longest to settle: ng 35 passes,
production 36 iterations.** The relation is asserted on **every cohort fixture** rather than on that
one, so it is an invariant of the port rather than an observation about a locus — and a real
difference in stopping point would break it by more than one. It is reported rather than asserted
at one sample, where production is also stopping on a different quantity and tying the two counts
together would pin a coincidence.

## 4. The fixtures, and what each can fail on

Reads that decide a genotype by ten nats cannot test a prior, so only the first of these is that
shape and it is there to catch wiring:

| fixture | what it can fail on |
|---|---|
| three samples, unambiguous reads | a table read in the wrong order, a sample's row handed to its neighbour |
| ten samples, one thin, cohort agrees with its reads | ten samples, a moving frequency, more than one pass |
| **ten samples, the cohort overturns the thin one's reads** | **a loop that ignored its prior** — the tenth sample's own reads say `0/1` and both loops call `0/0` |
| **one sample** | production's *other* E-step — the record-static one, which no cohort fixture reaches |
| the inbred panel | the recorded departure, in both directions |
| one sample of eight inbred, the rest not | the final pass reading sample 0's coefficient for every sample |
| three alleles | the alternative concentration **split**, asserted by value |
| eight samples the reads leave open | the prior carrying a locus outright, over 35 passes, with the pass-count relation |
| the two loops' threshold and cap | two independent constant pairs drifting apart |

**The third is the one that makes the oracle worth running.** Its assertion is that the call is
*not* the reads' own answer, computed from the table by an argmax the fixture does itself — so a
loop that dropped its prior would pass the others and fail that one.

## 4a. What the review found, and it is the eighth consecutive test that could not fail

**24 deliberate defects, 18 caught, 6 survived** — and the survivors named one accident: every
fixture ran at a diversity of 1 in 1,000 with read margins of six nats and up, so **nothing in the
*seed* ever decided a genotype** and only the leave-one-out cohort term did. Two of the six were
seed-only defects, and one of them was the division this file's three-allele fixture claimed to
test: deleting it changed no call. All six are caught by sibling tests elsewhere in ng's suite, so
the suite had no hole — this file had a reach gap and a false claim about itself.

**What changed.** The three-allele fixture now asserts the concentration **by value**,
`[1, θ/2, θ/2]`, which fails the moment the division goes; its doc says plainly that the genotype
comparison beside it cannot see the split and why. A **one-sample fixture** was added, because the
argument for excluding one sample was about the stopping rule and said nothing about genotypes —
and at one sample production runs a *different E-step*, which the oracle had never executed. A
fixture with **per-sample inbreeding coefficients** was added, since every other fixture gives every
sample the same one and so cannot tell this row's coefficient from row 0's; measured, it catches
the final pass making that substitution and **not** the frequency loop making it, which is stated
where somebody would otherwise rediscover it. The **threshold and cap** constants are now asserted
equal on the two sides, since they are two independent pairs that agree by coincidence of two
edits. And the module now says which three things it does *not* hold identical — the two constant
pairs and the transcendental backend — where it had claimed to hold everything.

**Five of 44 written claims were wrong**, every one a mechanism or a location rather than a number:
the "everything the loop reads" claim; "every fixture asserts" the pass relation, where the inbred
fixture asserted nothing about it; the three-allele split being "load-bearing"; a "thin sample the
prior has to place" whose own reads placed it; and a percentage quoted in the wrong direction.

## 5. Validation

- `cargo test --lib` → **4,923 passed / 0 failed / 14 ignored** (from 4,914 at the E2b commit; the
  nine are this module's).
- `cargo test --lib ng::calling::loop_parity` → 9 passed.
- `cargo fmt --all -- --check` and `cargo clippy --all-targets --all-features -- -D warnings` →
  both exit 0.

## 6. Banked for the owner

- **The oracle is fixed at more than one sample, and that is a real limit rather than a
  convenience.** At one sample production tests `max|p̂ − p̂'|` where ng tests copies over
  chromosomes; production's own comment explains why — a single-sample E-step is `p̂`-independent
  and reaches its fixed point in one iteration, so there is no trajectory to stop early on. The two
  rules therefore agree there for a reason, and comparing would compare two loops that have both
  already finished. **What is not covered is a difference that only shows at one sample**, and
  nothing here would see one.
