# ng calling loop — A1: the three shared types the loop takes and gives back

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step A1
**Design authority:** [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md),
[arch/calling_em_loop.md](../../ng/arch/calling_em_loop.md) §2
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Four category reviews raised one Blocker
> and thirteen Majors on the code below, and the fixes changed some of what §2 and §3 describe —
> `CallingScratch`'s type parameter lost its default, `prepare_for` became `prepare_for_locus`
> and takes the allele table, several names changed, and fifteen further tests landed. What was
> found, and what was done about it:
> [the review](../reviews/ng_calling_loop_a1_2026-08-25.md) and
> [the fixes](../reviews/fixes_applied_2026-08-25.md).

---

## 1. Plan

Add to [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs) the three types the
candidate-alleles and read-likelihood plans both deferred, now that every field they borrow
exists:

- **`CallingScratch`** — every buffer one locus's calling fills, allocated once per worker
  and reused at every locus, with candidate selection's fold buffers as a field.
- **`LocusEvidence`** — one locus's reads per sample, the one place the SNP/indel and
  repeat-tract paths meet.
- **`FrozenParameters`** — everything the parameter pre-pass fitted, gathered so one borrow
  crosses the calling seam.

Plus the two things the plan and the spec put in this file and nothing else could hold: the
**ordering contract**, asserted rather than promised, and a **carrier for the missing
genotype**.

Types only — no pass, no loop, no arithmetic.

## 2. Assumptions and deviations

Five, each a place where the shipped code was not the shape the architecture sketched. None
changes a design decision; each is recorded here because the architecture's sketch and the
code now differ.

### 2.1 `SampleGenotypeCall` becomes an enum

The architecture sketches `struct SampleGenotypeCall { genotype, genotype_quality }`
(arch §2), which cannot express a missing genotype — and
[`Genotype::new`](../../../../src/ng/types.rs) refuses an empty multiset by design. Spec §5.0
records that the ruling *"has a producer and no carrier"* and that whoever builds the loop
adds the variant here.

It is an enum rather than a struct with an optional quality because **a missing genotype has
no quality either**: the sample was never scored, so there is nothing for a quality to be the
confidence of. Emission must not write a missing `GT` with a poor `GQ` beside it, and an enum
is what makes that unwriteable.

**Built at A1 rather than at the plan's E3** — the end-to-end fixture where ng first calls
genotypes — which lists it among the three things that step has to join up correctly.
It is a shared type in the file A1 owns, and step C3 — the final pass, which fills
`LocusInference::per_sample` — cannot be written against the struct and then re-shaped later
without changing what C3 produces. Building it with the other shared types is one edit to one
file instead of two.

### 2.2 The missing-genotype flag rides on `LocusEvidence`, not on `GenericSampleEvidence`

[arch/read_likelihoods.md](../../ng/arch/read_likelihoods.md) §2.1 puts
`genotype_must_be_missing` on `GenericSampleEvidence`. **The shipped type does not carry it**
— its three fields are `supported`, `unmatched_q_sum` and `partials` — and
`src/ng/calling/likelihood/` is another branch's while this one runs.

It did not need to move there. Spec §5.0 rules that such a sample **leaves the loop before
the first pass** and is scored against no genotype at all, so the read likelihood never sees
it: a field on the evidence view would be one the row builder is handed and never reads. It
is instead a field of `GenericLocusSample`, a small pair of *(this sample's evidence, the
candidate step's ruling on it)* that `LocusEvidence::Generic` holds one of per run sample.

The pairing is deliberate rather than two parallel slices. The evidence and the ruling are
both per sample and both in run order, and a positional join between two separately-carried
lists is exactly the failure the spec's §5.0 closing paragraph names.

**`LocusEvidence::Ssr` carries no such field**, structurally — a repeat tract sets no sample
aside (spec §5.0.1), and there is nothing for the STR path to half-honour.

### 2.3 `CallingScratch` is generic over the repeat-tract emission model's scratch

The architecture writes `row_scratch: RowScratch` (arch §2), one thing. The likelihood
shipped two: `GenericRowScratch`, which is concrete, and `SsrRowScratch<ModelScratch>`, whose
parameter is the associated `Scratch` of whichever `SsrEmissionModel` scores tracts.

So `CallingScratch<SsrEmissionScratch = StutterSubstitutionScratch>` carries both, with a
default type so a run using the shipped model names no parameter. Guessing a concrete type
here would have picked one arm of a seam that exists to be swapped.

### 2.4 The prior's real API needs three per-allele buffers, not one

Architecture §2 lists one `concentration: Vec<f64>`. The genotype prior as built needs three
buffers of allele length at once, and they cannot be the same one:

- the **locus's seed**, filled once per locus by `fill_locus_concentration`, which returns a
  `Concentration<'_>` borrowing it;
- the **sample's own concentration**, filled per sample per pass by
  `fill_sample_concentration`, which reads the seed while writing this;
- the prior's **per-allele working space**, which `PriorRow::new` takes mutably while the
  concentration is borrowed immutably.

`CallingScratch` therefore has `seed_concentration`, `sample_concentration` and
`prior_allele_scratch`. It also has `prior_row: Vec<LogProb>` beside
`posterior_row: Vec<f64>`, because `PriorRow` writes `LogProb` and the posterior is a plain
probability.

### 2.5 `seed` is held by value; `ssr_strata` is never optional

`SpectrumSeed` is `Copy` and three `f64`-sized fields, so `FrozenParameters` holds it by
value where arch §2 writes `&'a SpectrumSeed`. `StratumFits` stays a borrow and is **not**
an `Option`: a run with no repeat tracts supplies `StratumFits::over(&[], BTreeMap::new())`,
whose lookups answer *no such stratum* — which is the honest answer, and one the STR row can
already act on.

## 3. Changes made

All in [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs).

| what | shape |
|---|---|
| `GenericLocusSample<'a>` | `{ evidence: GenericSampleEvidence<'a>, genotype_must_be_missing: bool }` |
| `LocusEvidence<'a>` | `Generic { region, per_sample: &[GenericLocusSample] }` / `Ssr { region, per_sample: &[SsrSampleEvidence], detail: &SsrDetail }`, both built through a constructor that refuses an empty sample list |
| `FrozenParameters<'a>` | private fields + `new` + accessors: calibration and contamination per read group, inbreeding per sample in run order, seed by value, `&StratumFits`, ploidy |
| `CallingScratch<SsrEmissionScratch>` | private buffers + `prepare_for` + accessors; holds `SelectionScratch`, `GenericRowScratch`, `SsrRowScratch<_>` |
| `SampleGenotypeCall` | struct → enum `Called { genotype, genotype_quality }` / `Missing`, with `genotype()`, `genotype_quality()`, `is_missing()` |

Three decisions inside those that are not restatements of the architecture:

**The ordering contract is a method.** `LocusEvidence::assert_agrees_with(alleles,
parameters)` checks two things in release: that the evidence's discriminant and the allele
table's `LocusKind` agree, and that the evidence's sample count equals the run's. Both are
caller bugs whose symptom is a wrong genotype rather than a crash — a repeat tract scored by
the SNP/indel row gives a plausible genotype at every sample, and two per-sample lists of
different lengths are two different orders.

**`prepare_for` takes the genotype table's view, not two integers.** At a diploid biallelic
locus the allele count is 2 and the genotype count 3; handed as bare integers, swapping them
leaves every buffer a legal length. Taking `&GenotypeTableView<'_>` means both come from one
object and cannot be crossed.

**Every buffer is poisoned with `NaN` on `prepare_for`, not merely resized.** `Vec::resize`
leaves the leading entries as they were, so a locus of the same shape as the last one would
silently reuse the last one's likelihoods and priors. The buffers are cleared and refilled.

**Flat tables are read through one indexer.** `lg_table` is `samples × genotypes` and
`per_sample_copies` is `samples × alleles`, both sample-major and both private;
`lg_row(sample)` and `sample_copies(sample)` are the only spellings of the index, and a
sample past the prepared count is refused by name.

## 4. Tests added

Fourteen, all in `src/ng/calling/mod.rs`'s test module. What each pins:

| test | what a wrong implementation would do |
|---|---|
| `frozen_parameters_refuse_a_contamination_list_of_another_read_group_count` | a mismatched pair is found at whichever locus first carries a read from the group past the end — or never, and then every genotype of the run is scored under somebody else's contamination |
| `frozen_parameters_take_an_empty_contamination_list_as_none_fitted` | absent contamination read as a fitted zero |
| `frozen_parameters_refuse_a_run_with_no_samples` | a run whose sample order went missing produces loci with no calls |
| `evidence_on_the_wrong_path_for_its_allele_table_is_refused` | a tract scored by the SNP/indel row: a different likelihood at every sample, and a plausible genotype at the end |
| `evidence_covering_a_different_sample_count_from_the_run_is_refused` | one sample's reads paired with another's inbreeding coefficient |
| `evidence_that_matches_its_locus_and_its_run_is_accepted_on_both_paths` | the accepting cell for a SNP/indel locus and for a repeat **tract**, and that the repeat variant carries no missing-genotype flag. **Not the repeat *bundle*, which is a third `LocusKind` and which this row originally read as covered** — the review's Blocker, now covered by two tests of its own |
| `evidence_naming_no_sample_at_all_is_refused` | evidence lost on the way in read as a locus nobody covered |
| `each_sample_reads_and_writes_its_own_row_of_the_scratch_tables` | a row sliced with the allele count instead of the genotype count reads a window straddling two samples' rows — every entry a real log-likelihood, none of them this sample's |
| `a_sample_past_the_prepared_count_is_refused` | reading the next sample's row, or a panic that names no sample |
| `preparing_a_locus_overwrites_the_previous_locus_of_the_same_shape` | the `Vec::resize` failure: locus *n* scored against locus *n−1*'s numbers |
| `advancing_makes_this_passs_copies_the_previous_passs` | the convergence test comparing a buffer against itself |
| `a_scratch_prepared_for_no_samples_is_refused` | as the parameters' own check, from the other side |
| `a_missing_call_carries_neither_a_genotype_nor_a_quality` | a quality beside a missing genotype, which conflates *not scored* with *scored and weak* |
| `a_locus_carries_a_missing_call_beside_called_ones_in_run_order` | the shape a cohort whose cap cut one sample's earned allele actually produces |

**The fixture that carries the most weight is the scratch's.** It writes each sample's own
index into that sample's row of both flat tables, so a wrong window names the sample it
actually came from rather than merely failing.

## 5. Validation

All in the container, from this worktree's own `scripts/dev.sh`.

| command | result |
|---|---|
| `cargo fmt --check` | exit 0 |
| `cargo clippy --lib --all-features --tests -- -D warnings` | exit 0, no warnings |
| `cargo test --lib` | **4,502 passed, 0 failed, 14 ignored** |

The baseline at the branch point was 4,488 passed / 0 failed / 14 ignored, so the fourteen
tests above are the whole of the difference.

`--all-targets` is deliberately not the clippy scope: it is red on `main` from pre-existing
lints in six examples and benches (`ng_joint_fit_harness`, `ng_duplicated_class_harness`,
`dhat_psp_reader`, `dhat_psp_writer`, `profile_posterior_engine`, `psp_writer_perf`,
`ng_joint_fit_perf`), none of them under `src/`.

## 6. Trade-offs and follow-ups

- **`prepare_for` poisons with `NaN`, which makes a forgotten write loud but not local.** A
  buffer some pass failed to fill fails at whichever check it is next handed to —
  `ExpectedAlleleCopies::new` refuses a `NaN`, `Concentration` checks its entries in debug —
  rather than at the omission. The alternative, zero-filling, would be a plausible value that
  reaches a genotype, which is the failure this file refuses everywhere else.
- **Not built here, and named so the next step does not have to find them:** the two
  `batches × alleles` contamination buffers, which the plan puts at E2a, where the contaminant
  frequency is wired per locus and per sample; the error-spread table the SNP/indel row takes,
  which arrives with B1, the first E-step; and `GenericEvidenceBuffer`, which the input edge at
  E1 needs in a concatenated form this type does not yet have.
- **One doc comment in another branch's module now reads as forward-looking when it is
  satisfied.** `SelectionScratch`'s comment in
  [src/ng/calling/allele_candidates/mod.rs](../../../../src/ng/calling/allele_candidates/mod.rs)
  says it *"becomes a field of `CallingScratch` when that type exists"*. It now is one. Left
  untouched: that module belongs to the candidate-alleles branch while it runs.
- **`arch/read_likelihoods.md` §2.1 still sketches `genotype_must_be_missing` on
  `GenericSampleEvidence`.** §2.2 above says why the loop does not need it there. Whether the
  architecture is amended is the owner's call, not this step's.
