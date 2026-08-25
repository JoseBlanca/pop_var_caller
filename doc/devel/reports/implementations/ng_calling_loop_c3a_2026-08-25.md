# ng calling loop — C3a: the two qualities the loop must compute while it still can

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step C3 — **the first of two commits**
**Design authority:** [spec/calling_quality.md](../../ng/spec/calling_quality.md) §3.1, §3.2, §4,
§5.1–§5.4, §7, §9, §10; [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md) §8
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Four agents returned **one Blocker that was
> a live wrong answer**, found independently by three of them with three different oracles: the
> fold read `ploidy` entries past what it zeroed, into the scratch's `NaN` fill, so the cohort's
> allele-count distribution was truncated to `0..=ploidy` **at every cohort of two or more
> samples**. At 63 diploid samples the site quality came out 46.3 Phred against production's
> 733.7 on the same table. The fix is one line. **Three of this report's own measurements were
> measurements of that defect**, and two of the module's tests were passing because of it.
> [The review](../reviews/ng_calling_loop_c3a_2026-08-25.md), and
> [what was done about it](../reviews/fixes_applied_2026-08-25_c3a.md).

---

## 1. Why C3 is two commits

The plan's C3 is one line — *"score every sample once more, take the highest-posterior genotype
and its confidence"* — and it is not one step. `spec/calling_quality.md` §3 puts three things
inside it, for one reason: **their inputs stop existing when the locus is released.** The
posterior row is a single reused buffer; the `samples × genotypes` likelihood table is per-worker
scratch; the per-sample read counts are borrowed from the merge.

So C3 splits, and the split is recorded here rather than made silently:

- **C3a, this commit** — the arithmetic, in a new module of its own: the per-sample genotype
  quality, the site quality before its artifact correction, and the type the nine pooled read
  counts will travel in.
- **C3b, next** — the final pass itself: score every sample once more, mint the owned genotype,
  gather the nine counts, and fill the two new `LocusInference` fields.

The split is along a seam that already existed. C3a is callable and testable standing alone;
C3b is wiring plus one open question (below).

## 2. What landed

**`src/ng/calling/quality/mod.rs`**, new:

- **`score_best_genotype(posterior_row) -> (GenotypeIdx, Phred)`** — the winner and
  `min(99, −10·log₁₀(1 − p_best))`, with the one-unit-in-the-last-place clamp and ties going to
  the lower index. One walk of the row returns both, because the caller needs both and the row
  is about to be overwritten.
- **`score_uncorrected_site_quality(...) -> Phred`** — *given every sample's reads, how unlikely
  is it that the cohort carries no copy of any non-reference allele?* Collapse each sample's row
  over genotypes into a row over non-reference copy counts; fold the samples into a distribution
  over the cohort's total, in the linear domain with per-sample rescaling; apply a Beta-Binomial
  prior; normalise and read the entry at zero.
- **`ArtifactTestCounts`** — the nine scalars §3.3 defines. A type only; C3b fills it.
- **`MAX_GENOTYPE_QUALITY` (99)** and **`MAX_SITE_QUALITY` (9999)**, both inherited and both
  marked soft.

**Four buffers on `CallingScratch`** and a `site_quality_buffers_mut` bundle, because nothing
here may allocate: the per-sample copy-count log-likelihoods, the two count-axis buffers the fold
alternates between, and the log-domain result. Sized from the genotype table's own ploidy, so
they cannot disagree with the table the fold walks.

### 2.1 The prior is the run's own spectrum, and that is the one intended difference

The collapse to (reference, any non-reference) turns the per-allele Dirichlet into a Beta on the
non-reference frequency, which induces a Beta-Binomial on the cohort's count. Production's two
concentrations are GATK constants, each carrying *"Revisit against the cohort calibration set"*
in its doc comment; nobody did. ng already holds the same two numbers fitted from the run's own
cohort, so it uses those (§5.4).

**Checked against §5.4's own table by a reviewer with an independent closed form: all fifteen
cells reproduce.**

### 2.2 Reuse rather than a third copy

The spec's reuse map says ng should write one `ln Γ` where production carries two. It needs
neither: `crate::genetics::lgamma` already exists and ng's genotype prior already uses it. Same
for the log-sum-of-two-exponentials the collapse needs — `log_sum_exp_2` was `pub(super)` inside
`genotype_prior` and is now `pub(crate)`, with its doc comment naming the third caller. Its
`−∞` guards are exactly what the collapse wants: a copy count no genotype reached stays
impossible rather than becoming a `NaN`.

## 3. Deviations

- **The step is split into two commits** (§1), recorded rather than silent.
- **`ArtifactTestCounts` keeps the spec's name**, though a reviewer is right that its first field
  is an `AlleleId` and not a count. Renaming a type the design document names is a bigger move
  than a naming fix; raised, not taken.
- **No new `Phred` constructor.** `src/ng/types.rs` invites one at "the step that first fills a
  `GQ` column", which is this one. It turned out not to be needed: both qualities are clamped
  into a known-good range before construction, so `try_new` plus a message is honest and adds no
  type surface.

## 4. What the review found, and what it cost

**The Blocker.** `fold_samples_into_allele_counts` zeroed `next[..ploidy + live]` and wrote
exactly that window — but the *following* sample, reading that buffer, tapped `ploidy` entries
further up, into the counts that sample had just made reachable. Those slots were never written,
and `prepare_for_locus` hands every scratch buffer over filled with `f64::NAN`. The `NaN`
multiplied in, survived the rescaling (`f64::max` returns the other operand), and was written out
as `−∞`.

**So from count `ploidy + 1` upward the cohort's allele-count distribution was silently declared
impossible, at every cohort of two or more samples.** Measured, three ways independently:

| oracle | result |
|---|---|
| brute-force enumeration of every cohort-wide genotype assignment | 157 of 252 cases disagree; **every case with two or more samples** |
| production's own `compute_qual_via_exact_af` on the same table and constants | 1 sample agrees to `1.3e-6` Phred; **63 samples: 46.3 against 733.7** |
| an independently written exact log-domain fold | **200 samples: 60.96 against 3316.55** |

**The fix is `next.fill(0.0)` beside the existing `current.fill(0.0)`.** With it, all three
oracles agree: brute force 252 of 252; the log-domain fold to a worst `1.15e-4` Phred over ploidy
1/2/4, one to five alternatives and up to 200 samples; production to a worst `1.23e-5` Phred.

**The same defect is latent in production** (`posterior_engine.rs:3624`), hidden on a freshly
grown scratch because `Vec::resize` zeroes, and reachable on a reused one — measured drift up to
**0.59 Phred at 8 samples**. Production's only test of that path is single-sample, which is
precisely the case the defect cannot reach. **Production is frozen and this branch does not touch
it**; the finding is recorded for its owner.

### Three of this report's measurements were measurements of the defect

| claimed | actually |
|---|---|
| 400 samples at 20 nats give 189.7 Phred, held down by the prior | 34,690 Phred uncorrected, **capped at 9999**; the prior costs 56 Phred, not 34,500 |
| the exact zero-term override "changes no answer on any fixture" | at 50 samples it is the difference between **4295.97 and the ceiling**; it starts mattering between 20 and 50 samples |
| 20 thin reference-looking samples move the shipped quality 0.032 → 0.237 | 0.032 → **6.55** |

**And the doc paragraph that recorded the second of these invited its own deletion** — *"a later
step trimming it should know that the tests will not object"* — of the line that is now holding
the answer together above 50 samples. That sentence is gone and a test stands where it was.

## 5. Tests

**Twenty-three.** Six on the genotype quality, thirteen on the site quality, four guards on the
buffer shapes.

The ones worth naming, all added or rewritten after the review:

| test | what it pins |
|---|---|
| `the_exact_zero_term_is_what_keeps_a_fifty_sample_cohort_off_the_ceiling` | 4295.97 with the override, 9999 without; unaffected at 20 samples |
| `the_collapse_strides_by_the_allele_count_and_not_by_the_ploidy` | a triallelic locus — **every other fixture is diploid biallelic, where the two strides coincide** |
| `shifting_every_likelihood_by_a_constant_leaves_the_quality_alone` | the running log scale, which every zero-peaked fixture leaves inert |
| `a_single_nan_beside_real_probabilities_is_refused` | the shape the winner alone cannot see — one `NaN` loses every comparison and the call comes out ordinary |
| `at_one_sample_the_site_quality_matches_the_closed_form` | an independently written Beta-Binomial oracle |
| `a_locus_nobody_carries_stays_at_nothing_however_many_samples_look_at_it` | 500 confident reference samples leave the quality at 1.3 in 100,000 of a Phred |

### ⚑ One claim of the specification did not reproduce, and it is the owner's

**`spec/calling_quality.md` §5.1 justifies the marginal formula by contrast**: the rejected
`Π_s P(hom-ref)` *"grows with cohort size at a site nobody carries"* where the marginal *"stays
bounded by what the few non-hom-ref samples actually justify"*.

**Measured on the fixed arithmetic, both grow, in the same proportion.** With deep reference
samples (12 nats), taking the cohort from 1 to 500 moves the rejected value from 0.0000267 to
0.0133 and the shipped one from 0.0000000267 to 0.0000134 — each 500 times its own base. With
thin ones (1 nat), the shipped quality grows *faster*: 0.0019 → 831.88 against 1.77 → 885.11.

The growth has a cause the section does not mention and that is not a defect: **in a cohort of
501 thin samples, "nobody carries this" is a far stronger claim than in a cohort of one, and thin
reads cannot support it.** The arithmetic is production's, unchanged, and agrees with it to
`1.2e-5` Phred — so this is a question about §5.1's *argument*, not about the code. The test that
was written against the contrast has been replaced by one asserting what is demonstrable.

## 6. Validation

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4643 passed; 0 failed; 14 ignored`. Before C3a: **4,620**.
- `cargo test --release --lib ng::calling --all-features` — `597 passed; 0 failed; 3 ignored`.
  Before C3a: **574**.
- **The release-held checks:** C3a adds nine. Downgraded all nine to `debug_assert` together and
  re-run under `--release`: `586 passed; 11 failed`, and **every one of the nine is reached** —
  the posterior row's emptiness, its total, the winning probability's range, the sample count,
  and five buffer shapes.

## 7. Follow-ups

- **C3b is the final pass**, and it carries one thing this commit could not settle: **the six
  strand and read-position counts the artifact summary needs are not on the calling evidence
  path.** `spec/calling_quality.md` §3.3 says they come *"from the merge's `SampleSupport`,
  borrowed through `LocusEvidence`"*. The merge's `AlleleSupport` does carry `num_fwd` and
  `placed_left` — but `GenericObservation`, the calling view of it, keeps only
  `(allele, read_group, num_reads, q_sum)`. Either the view widens by two `u32`s, or the counts
  are gathered at the input edge and passed in. **The view is where the spec puts them**, and
  widening it touches `src/ng/calling/likelihood/`, which is free now that its plan is merged.
  Raised here because it changes a type two plans consume.
- **`ArtifactTestCounts` is `pub` with nothing constructing it**, which is why it needed no
  dead-code expectation where every other new item did. C3b makes it real.
- **Production's `convolve_ac_linear` carries the same over-long read**, measured at up to 0.59
  Phred of drift on a reused scratch. Not this branch's to fix.
