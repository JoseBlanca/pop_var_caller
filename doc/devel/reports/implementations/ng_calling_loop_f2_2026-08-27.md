# ng calling loop — F2: the repeat-tract differential

**Step:** F2 of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — ng's loop at a repeat
tract under the existing caller's stopping rule and under its own.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §10.
**Date:** 2026-08-27. **Branch:** `ng-calling-loop`.

---

## 1. Why a differential and not an oracle

**The repeat-tract path does not get a parity oracle, and the reason §10 gives is about the two
loops rather than about either being wrong**: they converge on a different quantity at a different
scale, so two loops stopping at different points on one trajectory would disagree at any genotype
near a boundary, for no reason any document records.

*(An earlier draft of this report added that the existing caller's loop "cannot be handed ng's
table". It can:
[`em.rs`](../../../../src/ssr/cohort/em.rs)'s `run_pi_em` takes a per-sample × per-genotype
likelihood table, exactly as the SNP/indel oracle hands one to production's engine. What stops it is
module-private visibility, which the project's own rule would let a parity test widen. §10's reason
stands on its own and did not need the extra one.)*

**So §10 asks for something with a failing state instead:** run ng's loop under the existing
caller's convergence rule and tolerance, require the genotypes to match, then restore ng's rule and
report what moved.

## 2. The two stopping tests divide one movement by different totals

This is the step's one load-bearing reading, it was **wrong in the first draft**, and the way it
was wrong is worth keeping: the divisor was checked and the thing being divided was not.

Both loops take the largest per-allele change in the cohort's expected allele copies between
passes, and both turn it into a frequency before comparing it against a tolerance. **They divide by
different things:**

- the existing SSR caller adds its prior's pseudocounts to the copies *before* normalising —
  [`em.rs`](../../../../src/ssr/cohort/em.rs)'s `run_pi_em` opens each pass with
  `let mut expected = g0.to_vec()` — so its total is `chromosomes + pseudocount mass`;
- ng divides by the cohort's chromosomes alone.

**The first draft read `total = expected.iter().sum()`, concluded the divisor was the chromosome
count, and never checked where `expected` started.** It starts at the pseudocounts. So at one
nominal number ng's test is the **stricter** of the two, by a factor of one plus that mass over the
chromosomes: the SSR caller sees the same movement as smaller and stops sooner. Its own engine
records this as a real effect rather than a rounding one — a pseudocount-scaled readout *"does not
feed back"*, and testing it *"let a larger pseudocount damp the delta and stop the loop early"*
(spec §6).

**So what this step builds is ng's rule at two tolerances**, `1e-6` — the SSR caller's number —
and ng's own `1e-3`. **It is not a reproduction of that caller's rule**, which would need ng's
prior strength declared the counterpart of those pseudocounts: a claim about the two models that
nobody has made, and not one to make inside a test.

**The residual is measured, twice and independently.** Absorbing a plausible pseudocount mass into
the divisor — equivalently, loosening the tolerance by that factor — stops the tight arm at **four
passes rather than five**. One pass in five, on a comparison whose genotypes do not move either
way.

## 3. What it requires and what it reports

**Requires — and this is the failing state:** the genotypes match. On a tract of three samples at
**four reads apiece** the two tolerances give the same three calls, `1/1`, `0/1` and `0/0`.

**Four reads, with both neighbours measured before settling on it.** At twelve the loop settles in
one pass under the shipped tolerance and two under the tighter, so there is almost nothing to
compare; at two the *reads* stop deciding — the first sample comes out `0/1` rather than `1/1` — and
the fixture would be measuring the prior instead of the tolerance.

**Reports, asserted so the numbers cannot go stale:**

- **five passes against two.** The tighter rule costs two and a half times the work for the same
  three genotypes.
- **the two runs really did diverge.** Their cohort expected copies are not identical — asserted —
  which is what stops the genotype comparison being a comparison of one run with itself.
- **and the two stopping points land 4.8 × 10⁻⁵ of a chromosome apart** — asserted as a size
  within a factor of two, not as `> 0`, so that a last-place float wobble cannot pass for *the two
  runs diverged*. It is twenty times inside the looser tolerance. That is how far apart the two
  answers finished; it is **not** a promise either rule makes, since a convergence rule bounds the
  last step between two passes and says nothing about the distance to another rule's answer.

## 4. The second fixture, and a mechanism that had to be corrected twice

The same comparison at the tract ladder's **bottom rung** — a run whose repeat fit produced no
length spectrum anywhere, so every tract is seeded from a flat shape at one chromosome of belief.
Everything else is held at the first fixture's, so the one thing that changes is the prior's shape.
**Measured: three passes against one, where the fitted spectrum takes five against two.**

**Two wrong explanations preceded that number, and both were written before the measurement.** The
first said the loop would travel *further* from a weak prior; measured, it settles sooner. The
second said the strength of the prior was what drove the iterating; swept on this tract, raising a
fitted spectrum's concentration from 1 chromosome to 100 moves the pass count from 5 to 4 — barely,
and in the direction opposite to the claim.

**What separates the fixtures is the prior's *shape*.** An asymmetric spectrum pulls the
frequencies off the point the reads alone would put them, and the loop then has a trajectory to
walk; a flat one does not. Flattening the shape at a fixed concentration drops the count to 3 or 4.

**And the first version of this fixture reported one pass under both rules — the comparison
empty — because it changed two things at once**, the prior's shape *and* the substitution rates,
using a fixture helper whose own documentation says it is for a test that calls no tract. Holding
the rates gives the three against one above.

## 5. Validation

- `cargo test --lib` → **4,925 passed / 0 failed / 14 ignored** (from 4,923 at the F1 commit; the
  two are this step's).
- `cargo test --test ng_calling_loop_calls_genotypes` → 16 passed.
- `cargo fmt --all -- --check` and `cargo clippy --all-targets --all-features -- -D warnings` →
  both exit 0.

## 6. Banked for the owner

- **⚑ Whether ng's prior strength is the counterpart of the existing caller's pseudocounts.** If it
  is, this step could reproduce that caller's stopping rule exactly rather than run ng's at its
  tolerance — the difference is worth one pass in five here, and it grows as the cohort shrinks,
  since the pseudocount mass is a larger share of a smaller chromosome count. At one diploid sample
  it could be a factor of two. Nothing in the three calling documents identifies the two, and this
  step deliberately did not.

- **The differential is two tracts, both at three samples and four reads.** What it has not been run
  at is the other end of either committed axis — a thousand samples, or three hundred reads — where
  the frequency has less distance to travel and the two tolerances may stop in the same place for a
  different reason than the bottom rung's.
