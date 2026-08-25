# ng calling loop — C3b: the final pass, and the three things it is the last moment to compute

**Step:** C3b of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the second and last
commit of step C3.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §2, §5.0, §9;
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2;
[`spec/calling_quality.md`](../../ng/spec/calling_quality.md) §3, §4, §5, §6.3.
**Date:** 2026-08-25. **Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`.

---

## 1. What C3b is

C3a built the quality arithmetic in a module of its own and called none of it. **C3b is the pass
that calls it**, and the pass is the last moment three of the locus's outputs can be computed at
all:

- **each sample's genotype quality**, because `CallingScratch`'s posterior row is one
  genotype-length buffer that every sample in turn is scored into — when the loop returns, only
  the last sample's posterior still exists (`calling_quality.md` §3.1);
- **the site quality before its artifact correction**, because it folds the whole
  `samples × genotypes` likelihood table, which is per-worker scratch overwritten at the next
  locus (§3.2);
- **the nine pooled counts the correction consumes**, because eight of them are the evidence's,
  which is released with the locus, and the ninth needs the calls (§3.3).

Everything it produces goes into a `LocusInference`, including the cohort's expected copies
**as the loop left them** — recomputing those downstream from the called genotypes gives a
different number, because a call has already thrown away the uncertainty they carry
(`calling_em_loop.md` §9).

## 2. What landed

**One function and two helpers, all in `summarise_condition.rs`; two fields on
`LocusInference`; one buffer on `CallingScratch`.**

| item | what it is |
|---|---|
| `summarise_final_pass` | the pass: score every sample once more against the settled frequencies, take its winner and its quality in the same walk, mint the owned `Genotype`, pool the artifact counts, fold the site quality, assemble the `LocusInference` |
| `pool_reads_and_pick_primary_alternative` | pools every called sample's reads onto the alleles and picks the one the artifact tests weigh against the reference — the non-reference allele the most reads reached, `None` where there is no such allele. It is also the pass's one input-validation gate |
| `PooledArtifactCounts` | eight of the nine numbers while they are being summed: seven as `u64`, the eighth as an `f64` in the run's fixed sample order. The ninth is the allele the others are pooled *for*, and is never summed |
| `mint_genotype` | the table's copy-count row as an owned multiset — one allele per copy of the genome |
| `LocusInference::site_quality` | private, `pub(crate)` reader, one field written twice (§3.5) |
| `LocusInference::artifact_test_counts` | private, `Option<ArtifactTestCounts>`, public reader |
| `CallingScratch::pooled_allele_reads` | one `u64` per allele, so choosing the primary alternative costs no allocation |

`ArtifactTestCounts` had been `pub` with nothing constructing it since C3a; it is now real, and its
`expect(dead_code)`-free status is earned rather than an oversight.

### 2.1 Three decisions this step had to take, none of which the design fixes

**The artifact summary is an `Option`, and a quarter of loci are the reason.** `ArtifactTestCounts`
names a `primary_alternative: AlleleId`, and there are two ordinary situations with no such allele
to name: a locus whose candidate table came back as **the reference alone** — 27.4% of built loci
on the 63-accession tomato panel and 27.3% on HG002 at 30×, measured by candidate selection and
recorded on `SelectionVerdict::Selected` — and a locus whose alternatives no read reached. The
alternatives were to invent an id or to refuse the locus; production does neither, returning its
baseline unchanged in exactly these two cases ([`qual_refine.rs:79`](../../../src/vcf/qual_refine.rs)),
and `None` is that answer as a type rather than as an early return.

**§3.5's "one quality field" is enforced by visibility, and the assertion it asks for is the
ceiling check.** The rule has two halves: there is one quality field, and *nothing between the
worker and the correction stage may read it*. The first half is structural — one private field, one
reader, no second field. The second is enforced by that reader being `pub(crate)`: no consumer
outside the crate can see the uncorrected number at all, and the public accessor arrives with the
stage that makes the value public-worthy. What is left for an assertion is that the value came from
the function that owns the ceiling, which `LocusInference::new` now checks
(`site_quality <= MAX_SITE_QUALITY`). **A stronger form is available and was not built**: a state
marker saying whether the correction has run, which would let the emission gate assert it is
reading a corrected value — the exact defect §3.5 records. It is recommended to step 11's plan in
§7 rather than shipped here, because it would arrive with one of its two states unreachable.

**A sample's depth is its reads on the locus's alleles, and nothing else.** Two kinds of read are
outside it: a *partial* read showed no allele, so it stands behind none of the counts; and a read
whose allele candidate selection dropped reaches the calling view as pooled error mass with no
count beside it (`GenericSampleEvidence::unmatched_q_sum` carries a `q_sum` and no `num_reads`).
Production's depth is the same quantity — a sum over the alleles its record carries
([`qual_refine.rs:92`](../../../src/vcf/qual_refine.rs)).

### 2.2 What the pass deliberately does not do

**It does not re-sum the cohort.** It overwrites each *sample's* own expected copies, which is the
E-step's fourth step, and leaves the cohort row the loop converged on untouched. So the genotypes
are one E-step further on than the frequencies reported beside them — which is what production's
own final pass does, and at convergence the two differ by less than the threshold that stopped the
loop.

**It does not exclude a set-aside sample from the site quality.** It excludes such a sample from
the calls (`SampleGenotypeCall::Missing`, no quality beside it) and from the artifact counts and
their expectation, which is §6.3's rule. The site quality's fold runs over the scratch's whole
likelihood table, so its cohort is whichever samples the scratch was prepared for — and **whether a
set-aside sample gets a row there at all is D1's open choice**, the same one
`sum_cohort_expected_copies`' own note carries. Nothing wrong reaches a run today: such a sample's
likelihood row is never written, so the **loop's prior-free first pass** refuses the locus loudly
on the scratch's `NaN` sentinel — the E-step finds no usable score in a row of `NaN`s, before any
M-step has run — and the only way to reach the mismatch is a test that fills the row by hand.

## 3. Deviations from the plan

- **`LocusInference::new` grew from eight arguments to ten**, and every construction site was
  updated. The alternative — a second constructor, or a builder — would let a call site forget the
  two new fields, which is what the flat list prevents.
- **A fourteenth buffer on `CallingScratch`** (`pooled_allele_reads`), which arch §2's sketch does not
  list. It follows C3a's precedent, which added four for the site quality's fold, and the module's
  rule that a locus costs no allocation of its own. It is the only buffer there with no
  `UNWRITTEN_SCRATCH_VALUE` sentinel: reads are whole, so it is a `u64` in which a `NaN` is not
  expressible, and what replaces the sentinel is that its only reader zeroes and fills it in the
  same call.
- **`summarise_final_pass` takes nine arguments** and carries an `allow(clippy::too_many_arguments)`
  with its reason. Four of them are the seam's own — `LocusGenotyper::call_locus` takes the
  evidence, the parameters, the candidates and the scratch — and the five that remain are the prior
  model, the loop's outcome and the two warrants a locus carries, all of which D1 has in hand at the
  call site. All nine have distinct types, so a transposition is a compile error rather than a wrong
  answer.

## 4. Tests

**Twenty-eight**, measured as the difference the commit makes to the library target:
4,644 → 4,672. Twenty-two in `summarise_condition.rs`, six in `calling/mod.rs`. **Six of the
twenty-eight are the review's**, and every one of those six closes a test that could not fail
(§7).

**What the fixtures were made to disagree about.** C3a's review found three mutations surviving a
green suite, all hidden by one habit: every fixture in the module was diploid, biallelic, and had
rows peaking at zero. The fixtures here vary five things on purpose, and each variation has a test
that fails without it:

| the shape varied | the test | what a fixture without it hides |
|---|---|---|
| **three alleles, and the winner is the second alternative** | `the_primary_alternative_is_the_allele_the_most_reads_reached_over_two_read_groups` | a hard-coded `AlleleId(1)` passes every biallelic fixture |
| **two read groups per allele** | the same test | a fold that reads only the first row of a sample |
| **ploidy four** | `at_ploidy_four_a_call_carrying_one_copy_expects_a_quarter_of_its_reads` | a divisor of two, and a mint that assumes two copies — at ploidy 2 the ploidy, the copies of a heterozygote and the alternative count are all 2 |
| **one allele, the reference alone** | `a_locus_called_over_the_reference_alone_carries_no_artifact_summary` | the `None` arm, on a quarter of built loci |
| **a sample set aside** | `a_sample_the_candidate_step_set_aside_is_missing_and_counted_nowhere` | counting an uncallable sample's reads: its 20 alternative reads would take the pool from 6 of 12 to 26 of 38 |

The rest, and what each pins:

| test | what it pins |
|---|---|
| `one_samples_call_and_its_confidence_match_the_arithmetic_done_by_hand` | posterior `1/7, 2/7, 4/7` → genotype `1/1` at `10·log₁₀(7/3)` = **3.6797678**, and all nine counts by hand |
| `a_tie_between_two_alternatives_goes_to_the_lower_allele_id` | the first *strict* maximum, so a run reproduces itself |
| `the_expected_alternative_reads_come_from_the_calls_and_not_from_the_frequency` | the calls give exactly 4; the fitted frequency (0.23898 of four chromosomes) would give 3.8237 |
| `the_cohort_expected_copies_are_the_loops_and_not_recomputed_from_the_calls` | `[2.6759, 1.3241]` on the record against the `[2, 2]` two heterozygous calls imply — two thirds of a copy apart, which is twice the third of its probability each sample keeps on the homozygous reference |
| `the_site_quality_on_the_record_is_the_fold_of_the_table_the_loop_built` | 42.373 Phred, against a second independently prepared scratch |
| `reads_that_saw_only_part_of_the_locus_are_in_neither_the_depth_nor_the_total` | a 50-read partial beside the 8 counted reads — a factor, not a rounding |
| `a_repeat_tract_carries_a_site_quality_and_no_artifact_summary` | §8's ruling: the strand and position tests are the SNP/indel path's |
| `a_locus_whose_alternatives_no_read_reached_carries_no_artifact_summary` | the other `None` arm; the quality is still computed, at 8 in 100,000 of a Phred |

## 5. Validation

All in the container, from this worktree.

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4672 passed; 0 failed; 14 ignored`. Before C3b: **4,644**.
- `cargo test --release --lib ng::calling --all-features` — `626 passed; 0 failed; 3 ignored`.
  C3a recorded 597 and `5781a6a9` added one, so this is the same **+28**.
- **The release-held checks: C3b adds seven.** Downgraded all seven to `debug_assert` together and
  re-ran under `--release`: `618 passed; 8 failed`, and **every one of the seven is reached** — the
  evidence's sample count, the inbreeding coefficients' length, an observation's allele being in
  range, its strand and placed-left counts being shares of its own reads (two tests, one check),
  the site quality's ceiling, and the artifact summary's primary alternative being neither the
  reference nor an allele the locus lacks. The files were restored from byte-identical backups and
  the release suite re-run green at 620 before anything was committed.

## 6. What the review found, and what it cost

**Three agents, each in its own worktree: one on tests and mutation, one on naming/errors/idiomatic
/smells/refactor-safety, one briefed only to re-derive the diff's own numbers.** Verdict:
**2 Blockers, 4 Majors, 5 Minors, 1 Nit**, plus **10 wrong claims of 42 checked**. Every finding
was applied.

**Both Blockers were tests that could not fail, and both were the same fixture habit the plan
warned about.** No fixture asserting a genotype quality had more than one called sample — so
overwriting every sample's quality with the first sample's left the whole 4,666-test suite green,
which is exactly the defect `calling_quality.md` §3.1 says this pass exists to prevent. And the one
three-allele fixture, built to defeat a hard-coded `AlleleId(1)` in the *choice* of alternative,
asserted eight of the nine counts but not the ninth — so hard-coding `copies_of_each_allele[1]` in
the expectation also passed everything. Both are now pinned; the second cost one added assertion.

**Two of the four Majors were the same shape.** The reference row's forward and placed-left counts
were equal in both fixtures that asserted them, so swapping the two accumulations survived — the
hand-computed fixture's reference row is now 3 reads, 2 forward, 1 placed left. And the set-aside
sample's exclusion from the *choice* of primary alternative was untested, because that fixture was
biallelic where there is only one alternative to choose.

**Ten of forty-two quantitative or mechanism claims were wrong, and the three that mattered were
mechanisms rather than numbers**, which is the same pattern C3a's review found. This report said a
set-aside sample's unwritten likelihood row is refused by the M-step; it is refused by the loop's
prior-free *first pass*, in the E-step, before any M-step runs. A doc comment said each sample of
one fixture keeps a fifth of its probability on the homozygous reference where the posterior is
`[0.3438, 0.6507, 0.0055]` — a third, and the third is what makes the two-thirds-of-a-copy gap
arithmetic close. **And a doc comment repeated the plan's claim that nothing in ng can score a
repeat tract because the read-likelihood row is `unimplemented!()` — which is no longer true**: the
row shipped with that plan's H1 and H2, and the only `unimplemented!()` left under
`src/ng/calling/` is a `#[cfg(test)]` oracle's. What still blocks a tract end to end is the
repeat-tract *candidate* path, which is unwritten. The plan's blocker note is corrected in this
commit.

Everything in §5's validation counts was checked and correct, including the battery.

## 7. Follow-ups

- **D1 owns the set-aside sample's row**, and now owes it in two places rather than one:
  `sum_cohort_expected_copies` (recorded at B2) and the site quality's cohort (§2.2 above). Both
  are the same choice — skip those rows, or never give them one.
- **Step 11's plan should consider a correction-state marker on `LocusInference`.** Not a second
  quality field, which §3.5 forbids: a marker saying whether the field holds the baseline or the
  corrected value, so the emission gate can assert it is reading the corrected one. That is exactly
  the defect §3.5 records, and the marker is the only form of the invariant that an assertion can
  reach. It is not built here because one of its two states would be unreachable until that stage
  exists.
- **`weakest_provenance` and `seed_diversity_unreachable` arrive as arguments.** Nothing computes
  them yet; the row builders and the prior will, and D1 is where they are gathered.
