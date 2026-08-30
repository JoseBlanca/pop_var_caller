# ng call quality (step 13) — implementation plan

**Status:** draft, 2026-08-30. The build order for the **artifact correction** — the third of the
three numbers [`calling_quality.md`](../spec/calling_quality.md) defines, and the only one not yet
built — plus the production parity oracle the other two were shipped without. Design is settled in
that spec; there is no separate architecture document, because the spec's §10 carries the types and
§11 carries the reuse map. This plan turns that design into build order and is **not** a place for
new design: every question it met is resolved in the spec, and the two the spec leaves open (§13's
Q1 and Q2) are handled the way the spec says to handle them — Q1 by porting unchanged and measuring
later, Q2 by the coder at the keyboard against a measured agreement (step B2).

**Two of the three numbers are already built**, in step C3b of
[`calling_loop.md`](calling_loop.md) (2026-08-25):
[`quality/mod.rs`](../../../../src/ng/calling/quality/mod.rs) computes the per-sample genotype
quality and the site quality before correction, and
[`summarise_condition.rs`](../../../../src/ng/calling/inference/summarise_condition.rs) pools the
nine-number artifact summary in the worker. What is missing is the arithmetic that turns the
baseline into the quality a file carries, and the oracle that says the baseline is right.

---

## Scope

**In:** `src/ng/calling/quality/artifact_correction.rs` and `src/ng/calling/quality_parity.rs`;
`ArtifactPenalties` and the two ramp constants of spec §10; the regularised-incomplete-beta
two-sided binomial tail; the allele-balance penalty with its two guards; the strand and
read-position penalty with its power ramp; the corrected quality, floored at zero; and the
differential against production's `refine_qual` over the same nine numbers.

**Out (later plans, each with an owner):**

- **The stream wiring** — which stage runs the correction, the in-place overwrite of
  `LocusInference`'s one quality field, and the emission threshold that reads it. Spec §3.4's
  Decision hands the wiring to **step 11's document, unwritten**; this plan supplies the penalties
  and the corrected quality *as functions of the summary* and touches no stream code. The
  `pub(crate) uncorrected_site_quality()` accessor already carries a lint-allow naming that step.
- **The GIAB depth-ladder measurement** — spec §14's tests 6 and 7 (the false-positive quality must
  not rise with depth; low-depth recall must not fall) and §13's **Q1** (do the two tests hold on a
  cohort, and are the ramp endpoints right at three reads a position). Both need ng to call a real
  cohort and emit a file, which nothing does yet. **Home:** step 11's plan, or its own measurement
  note once a run exists. Until then the ramp endpoints ship at production's `(3, 7)` with their
  provenance in the doc comment, which is the spec's stated leaning.
- **Worker-count bitwise invariance** — spec §14 test 5's first half. The fold this plan touches
  runs inside one worker and sees no thread count; the test belongs where the run does. **Home:**
  step 11's plan.
- **Bias annotations as VCF fields in their own right** — spec §12. **Home:** the VCF output
  document.
- **Repeat tracts** — spec §8 says nothing here is written for them. **Home:** the repeat-tract
  sibling of that section.
- **A genotype quality that knows the site is an artifact** — spec §1.2 excludes it; it would
  change calls, not just their confidence. **Home:** its own investigation, if a measurement asks.

---

## Principles (how the order was chosen)

- **The oracle before the code it judges.** Milestone A builds the production parity for the site
  quality baseline **first**, even though that baseline shipped three days ago. The correction is a
  *subtraction from* that baseline, so its own differential (Milestone D) compares two numbers that
  each rest on it; a baseline nobody has checked against production would let a disagreement in D
  be blamed on either half. This is spec §14's test 2, owed by C3b and unpaid.
- **Reuse over rewrite.** This is a port. It calls ng's existing `crate::genetics::lgamma` rather
  than transcribing production's hand-written Lanczos approximation (spec §11), and it
  reproduces production's guards, clamps and ramp rather than re-deriving them. No test's threshold
  is re-argued here — spec §6.2 records where each came from.
- **The algorithmic heart before the plumbing.** The binomial tail (B2) is built and tested against
  production's own discrete sum before either test that calls it. The two tests are then pure
  functions of nine scalars, and the assembly (C3) is a subtraction.
- **Isolate the step whose failure is silent.** The tail returns a plausible number for a wrong
  answer: a mis-transcribed continued fraction gives a Phred that is finite, positive and
  ordinary-looking, and it reaches the QUAL column as a slightly wrong confidence rather than as a
  panic. **B2 lands as its own commit, not bundled into B1 or C1**, with its oracle green before
  and after, so a bisect can find it if a quality moves later.
- **Types first, then implementation**, within every milestone (project rule).
- **Verify against ground truth.** The north-star test is agreement with production's
  `refine_qual` on the same nine numbers (D1), not self-consistency — a port is correct when it
  matches what it ports.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Container builds / ungated.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); `ng` stays a plain
  module behind no feature gate.

---

## Preconditions (already in place)

- **Step C3b of [`calling_loop.md`](calling_loop.md) is done**, so all three of the correction's
  inputs exist: `score_uncorrected_site_quality` writes the baseline,
  `summarise_condition::into_summary` fills every field of `ArtifactTestCounts` exactly as spec §10
  declares it, and `LocusInference` carries both — the quality behind
  `uncorrected_site_quality()`, the counts behind `artifact_test_counts()`, the second already
  returning `None` for the two cases production returns its baseline unchanged for.
- **`crate::genetics::lgamma`** — already used by `score_uncorrected_site_quality`, so the tail
  needs no second one. **It is `libm`'s, where production's is a hand-written Lanczos
  approximation**, so a port that reads it is not bit-identical to its source by construction; the
  two agree to about one part in `1e-13`. Spec §11 expected ng to write the Lanczos; it has
  something better, and D1 measures what the difference costs.
- **`Phred`** ([`types.rs:429`](../../../../src/ng/types.rs)) refuses negatives, `NaN` and
  infinities at its one constructor, and normalises `-0.0`. Every penalty and the corrected quality
  cross into it.
- **The production reuse target and oracle:**
  [`qual_refine.rs`](../../../../src/vcf/qual_refine.rs) (`refine_qual` and the four functions
  under it) and [`allele_balance.rs`](../../../../src/var_calling/allele_balance.rs).
- **Milestone A needs no change to production.** `run_em_columnar` is `pub(crate)` and the
  `EmOutputs` it returns carries `qual_phred`, which is `compute_qual_via_exact_af`'s return over
  the *input* likelihood table and the record-static pseudocounts
  ([`posterior_engine.rs:2479`](../../../../src/var_calling/posterior_engine.rs)) — the same two
  inputs ng's function takes, and neither of them moved by the EM.
  [`loop_parity.rs`](../../../../src/ng/calling/loop_parity.rs) already drives that entry point and
  its fixture shape is the one to copy.
- **⚠ The prior the site quality reads is not the prior the frequency loop reads**, and a fixture
  built on the second tests nothing. Production carries two concentration pairs per record:
  `scratch.alpha`, `[1, θ̂/k, …]` from the run's nucleotide diversity, which the EM uses and which
  `loop_parity` holds identical; and `scratch.pseudocounts`, `[10, 0.01, …]` from four compiled-in
  GATK constants, which the site quality uses and which spec §5.4 says ng replaces.
  `with_nucleotide_diversity` moves the first and not the second.
- **Milestone D needs one visibility widening, and it is the freeze exception** (CLAUDE.md; widen
  `pub(crate)`, change nothing else): `mod qual_refine` in
  [`vcf/mod.rs:38`](../../../../src/vcf/mod.rs) is private and `refine_qual` is `pub(super)`. Both
  become `pub(crate)`. No production logic, constant or signature moves. `PosteriorRecord` already
  implements `VcfWritable`, so the differential builds one of those rather than writing a fake with
  twenty stub methods.

---

## The steps

### Milestone A — the oracle for the baseline that already shipped

**A1. The site quality against production's, at production's constants.**  ✅
New `src/ng/calling/quality_parity.rs`, beside `loop_parity.rs` and `genotype_table_parity.rs` and
named for the same reason. One fixture's genotype log-likelihood table, in the genotype order
`genotype_table_parity` already pins, goes to both sides: to `run_em_columnar`, whose
`EmOutputs.qual_phred` is production's answer, and to ng's `score_uncorrected_site_quality`. **ng
is seeded from `DEFAULT_REF_PSEUDOCOUNT` and `DEFAULT_SNP_ALT_PSEUDOCOUNT` — the site quality's
pair, not the EM's — so the two priors are the same construction rather than two transcriptions**:
production sums the per-allele pseudocount over the alternatives to reach the Beta's `α_alt`, and a
hand-typed total would agree at two alleles, part at three, and go on passing after a change to
production's constant. One fixture moves both constants, since every other one runs at the shipped
values where a port that hard-coded them would pass. Fixtures at one, twelve and forty samples, one
triallelic, one with no evidence at all, and one driven past ng's ceiling — which production has
none of, so it is the one place the two are asserted to part. *Depends:* none. *Source:* spec §11,
§14 test 2.

**A2. The second arm, and the permutation tolerance.**  ✅ *(landed in A1's commit — see below.)*
Same fixtures with ng's **fitted** seed instead of production's constants: the movement is reported
and asserted in sign against spec §5.4's table, at both ends of it. **A silent agreement here is a
failure**, not a pass — it would mean the seed never reached the prior — so the test asserts a
non-zero difference in the stated direction rather than a bound. Alongside it, spec §14 test 5's
second half: permuting the cohort moves the fold's summation order, so the same locus under a
permuted sample order must agree to a tolerance, asserted as a tolerance and never as bitwise
equality. *Depends:* A1. *Source:* spec §5.4, §11, §14 tests 2 and 5.

**Not split from A1, and the reason is the module's own prose.** A commit holding arm one alone
would ship a file whose doc comment describes two arms and contains one, and whose central claim —
*only one of these two is parity* — would have nothing to point at. The two arms are seven tests
and two in one 470-line module; splitting them buys no bisectability, because neither arm can fail
silently: an arm-one break is a numeric disagreement and an arm-two break is a direction.

> **Checkpoint A:** ng's site quality reproduces production's number where the two priors are the
> same, and the fitted seed's departure from it is measured rather than assumed. Pause for review.

### Milestone B — the types, and the tail (the numerical heart)

**B1. The types and the two constants.**  ✅
In a new `src/ng/calling/quality/artifact_correction.rs`, declared from `quality/mod.rs`:
`ArtifactPenalties { allele_balance: Phred, strand_and_read_position: Phred }` and the two
`pub const`s `BIAS_RAMP_NO_POWER_BELOW = 3.0` / `BIAS_RAMP_FULL_POWER_AT = 7.0`, each with spec
§6.2's provenance in its doc comment — read off the GIAB HG002 alternative-read distributions at
one sample in June 2026, soft, and Q1's to revisit. **Typed constants and not production's
`PVC_BIAS_RAMP` environment variable**, which is the shape spec §3.5 exists to prevent. A new file
rather than more of `quality/mod.rs`, which is at 1,061 lines and owns a different input — a
`samples × genotypes` table, where everything here is nine scalars (module_layout §Organizing
principles, rule 3's "extract a coherent chunk when it tells you to"). No logic. *Depends:* none.
*Source:* spec §10, §6.2, §3.5.

**Three constants rather than two**, because production's allele-balance guard is an inline
`0.9` and the same rule applies to it: `ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE`. And the two range
properties are `const { assert! }` items rather than tests — they are properties of literals, so
the compiler settles them where they are written.

**B2. The two-sided binomial tail.**  ✅ *(its own commit, as required.)*
The regularised incomplete beta and the tail built on it, ported from
[`qual_refine.rs:305-457`](../../../../src/vcf/qual_refine.rs) (`betacf`, `reg_incomplete_beta`,
`binom_cdf_le`, `binom_sf_ge`, `binom_two_sided_p_beta`, `tail_phred`), reading ng's
`crate::genetics::lgamma` instead of that file's private `ln_gamma`. **This settles spec §13's Q2 at
its leaning: one implementation, the closed form.** Production keeps a discrete sum below 2,000
reads only to hold its own output byte-identical at depths it had already validated; ng has no such
obligation, its oracle is a differential, and a cohort's pooled read count is above 2,000 in any
case — so the sum would be dead code on the runs that matter. **The sum is ported into
`#[cfg(test)]` instead**, where it is this step's oracle: the two must agree across a grid of
`(k, n, p)` spanning both tails, both boundaries and the degenerate `n = 0`, and the step records
the tolerance they agree to. *Depends:* B1. *Source:* spec §11, §13 Q2.

**Measured: the two agree to `7.0e-13` across 8,155 comparisons** — read totals from 1 to 999,
expected shares from 1 in 100 to 99 in 100, every outcome from none to all. **And ng's `lgamma` is
`libm`'s where production's is a hand-written Lanczos approximation**, so the port is not
bit-identical to its source by construction; what that costs at the end of the whole correction is
D1's to measure.

> **Checkpoint B:** the tail agrees with an exact discrete sum to a stated tolerance, and the two
> ramp endpoints exist as typed constants carrying their provenance. Pause for review.

### Milestone C — the two tests, and the subtraction

**C1. The allele-balance penalty.**  ✅
A function of `ArtifactTestCounts`: the expected alternative-read fraction is
`genotype_expected_alternative_reads / total_reads`, clamped as production clamps it, and the
penalty is the tail of the observed split against it. **Both guards, and both are load-bearing
rather than defensive.** Only a *deficit* is charged — an excess is a different phenomenon this
test says nothing about — and the test is skipped entirely where the called genotypes expect a
fraction at or above 0.9, because a homozygous-variant sample's few reference reads are sequencing
error and a binomial against a probability near one reads them as a deficit. Unit tests: a
well-balanced heterozygote pays about nothing; a 20%-of-50% split at depth pays and pays more at
twice the depth; an excess pays zero; a cohort at 0.95 expected pays zero. *Depends:* B2.
*Source:* spec §6.2, §6.3.

**Measured, and the ratio is the point:** one read in five where the genotypes say half should
costs **46.2** Phred at 50 reads and **430.8** at 500. A penalty that did not grow with depth
would be swamped by a site quality that does. **And neither end of the expected-share clamp can
bind** while the deficit rule and the 0.9 guard stand — the first puts the share above
`1 / total_reads` wherever the tail is reached, the second returns before the ceiling. Both are
kept as production has them, as the net under two constants somebody may move.

**C2. The strand and read-position penalty, with its ramp.**  ✅
The larger of two tails — the alternative reads' forward-strand fraction, and their placed-left
fraction — each taken against **the reference reads' own fraction at the same site**, clamped to
`[0.01, 0.99]`, falling back to one half where there are no reference reads. Using the reference
reads as the expectation rather than a fixed half is what keeps the test honest at a site whose
coverage is one-sided for innocent reasons. The raw penalty is then scaled by the ramp of B1:
nothing at or below three alternative reads, full at seven or more, linear between. Unit tests:
an evenly-sampled alternative pays about nothing; a one-strand pile-up at 40 reads pays; **the same
pile-up at 3 reads pays exactly zero and at 5 reads pays half**, which is the property the ramp was
added for and the one that a transcription error in the ramp leaves silently missing. *Depends:*
B2. *Source:* spec §6.2.

**C3. The corrected quality.**  ✅
One function: baseline `Phred` and `ArtifactTestCounts` in, `(Phred, ArtifactPenalties)` out, the
quality being `baseline − allele_balance − strand_and_read_position` floored at zero. `None` counts
never reach it — `LocusInference` already carries `Option<ArtifactTestCounts>` and returns `None`
for the two cases production returns its baseline unchanged for — so this function takes the
summary by value and the *absent* case is the caller's, which is the stream stage's and not this
plan's. Unit tests: penalties exceeding the baseline floor at zero rather than going negative
through `Phred::try_new`; a clean locus comes back with its baseline and two zero penalties; the
uncorrected value is recoverable as the sum wherever the floor did not bind, which is spec §3.5's
claim and the reason there is no second quality field. *Depends:* C1, C2. *Source:* spec §2, §3.5,
§6.

> **Checkpoint C:** the correction is complete as arithmetic over nine numbers, each test's guards
> and ramp pinned by a fixture that fails without them. Pause for review.

### Milestone D — the differential against production

**D1. ng's correction against `refine_qual`, over the same nine numbers.**  ☐
In `quality_parity.rs`. Build a `PosteriorRecord` whose per-sample `scalars` and `best_genotype`
pool to a chosen `ArtifactTestCounts`, hand it and a baseline to production's `refine_qual`, hand
the same summary and the same baseline to C3, and require the two corrected qualities to agree.
The widening of the preconditions lands here, in its own commit, touching two lines of
`src/vcf/`. Cases: a clean biallelic heterozygote, a strand-piled artifact above and below the
ramp, an allele-balance deficit at two depths, a homozygous-variant site the 0.9 guard skips, and a
locus with no alternative reads at all — where production returns the baseline unchanged and ng's
caller passes `None`. **`PVC_BIAS_RAMP` must be unset for the run**, since production reads it once
into a `OnceLock` and a set value would silently make the two ramps different; the test asserts
that rather than assuming it. *Depends:* C3, A2. *Source:* spec §11, §14 test 2.

**D2. What the correction charges, on a real cohort's numbers.**  ☐
Not a benchmark and not a measurement of recall — those are out of scope above, and need a run that
does not exist. This is the smallest honest thing available now: over the artifact summaries a
tomato-panel and an HG002 fixture locus set already produce through the calling loop, report the
distribution of each penalty and how often each guard and the ramp bind. It answers one question
Q1 will ask later — whether a cohort's pooled alternative-read count crosses the ramp's seven for
reasons that have nothing to do with one sample's power — and it costs one example program.
Records the numbers; changes no constant. *Depends:* D1. *Source:* spec §13 Q1.

> **Checkpoint D:** ng's corrected quality reproduces production's on the same inputs, and what the
> two penalties charge on this repository's own cohorts is written down. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| **A** — the baseline's oracle | ng's `score_uncorrected_site_quality` against `PosteriorRecord.qual_phred` from `run_em_columnar` on the same likelihood table, ng seeded at `(ALPHA_REF, θ)` so the priors are one construction; second arm at the fitted seed asserted non-zero and in spec §5.4's direction; cohort permutation to `1e-6` |
| **B** — the tail | the closed-form incomplete-beta tail against a `#[cfg(test)]` port of production's exact discrete sum, over a `(k, n, p)` grid covering both tails, both boundaries and `n = 0`, to a recorded tolerance |
| **C** — the two tests | unit fixtures that fail without each guard: deficit-only, the 0.9 skip, the reference-fraction expectation, and the ramp charging zero at 3 alternative reads and half at 5 |
| **D** — the port | agreement with production's `refine_qual` on the same nine numbers across six locus shapes, with `PVC_BIAS_RAMP` asserted unset; plus a recorded distribution of both penalties over tomato and HG002 fixture loci |

---

## Out of scope (next plans)

- **Step 11 — the output stream.** The stage that runs C3, the in-place overwrite of
  `LocusInference`'s single quality field (spec §3.5), the emission threshold reading the field the
  stage wrote, and spec §14's tests 6, 7 and 8. Its document is unwritten; this plan's C3 is the
  one function that stage calls.
- **The depth-ladder measurement and §13's Q1.** The false-positive quality curve, low-depth
  recall, and whether `(3, 7)` holds on a 63-accession panel at three reads a position. Needs a run
  end to end; **home:** step 11's plan or a measurement note beside it. D2 supplies the part that
  does not need a run.
- **The repeat-tract sibling of spec §8**, which owns whatever a tract's quality turns out to be.
- **Bias annotations as VCF fields** (spec §12) — the VCF output document.
