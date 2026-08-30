# ng calling loop — E2a: the contaminant frequency, per locus and per sample

**Step:** E2a of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the second half of the
contamination mixture, and the batching it is drawn against.
**Design authority:** [`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.6 and
§6.1; [`arch/parameter_prepass_joint_fit.md`](../../ng/arch/parameter_prepass_joint_fit.md) §1.6;
[`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §2, §3, §8, §13.
**Owner's ruling this step implements:** the plan's E2a entry, 2026-08-26 — the genotype-likelihood
table's build splits in two.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

A run can now be called with a contamination fraction fitted. The driver's refusal is gone. Three
things had to exist for that: **who was sequenced beside whom**, which is declared by the user and
was specified but unbuilt; **the two fills** that turn the loop's own expected allele copies into
the frequency an observation's allele has in the contaminating population; and **a split of the
genotype-likelihood table's build**, because that frequency moves with the loop and the table the
loop reads no longer stands still.

## 2. The split, and what it keeps

**The contradiction it resolves.** `spec/read_likelihoods.md` §3.6 says that with contamination on
*"a caller may no longer cache a whole row across iterations"*. Step D2 asserts the opposite — the
table built once per locus. Both are right about different halves, and the owner's ruling separates
them:

- the **emission** — the answer to how one copy of one allele produced one observed sequence —
  reads no allele frequency, so it is still computed once per `(sample, observation, candidate)` per
  locus;
- the per-genotype **fold** over those emissions reads `q(o)`, so it runs again whenever the
  frequencies move, which is once a pass — **and only where a fraction was fitted**.

**What the emission actually is on the SNP/indel path**, since this is the thing that had to be
made cacheable and nothing had named it before: one *charged error* `ε̄` per observation — the
geometric mean of its reads' own error probabilities, scaled by the read group's calibration — and
one *compatibility verdict* per `(partial read, candidate allele)`, which is the byte comparison
deciding whether that allele could have produced a read that saw only part of the locus. Both are
now filled by `fill_generic_emissions` into a `GenericRowScratch` **per scratch row**, where before
there was one such scratch reused across samples. `assemble_genotype_log_likelihood_row` is the
fold. `genotype_log_likelihood_row` still exists with its old signature and is now the two called in
order, so every test written against the composed row is unchanged — and
`the_composed_row_is_its_two_halves_called_in_order` asserts the two spellings agree to the bit.

**What `EmissionCost` counts now says which half.** `table_builds` and `row_builds` are gone;
`emission_builds`, `emission_row_fills`, `emission_evaluations` count the expensive half and
`table_assemblies`, `row_assemblies` the cheap one. **D2's numbers did not move**: on an
uncontaminated run every fixture reports `emission_builds: 1` and `table_assemblies: 1` where it
reported `table_builds: 1`, with `emission_evaluations` unchanged, because with no frequency in the
formula the folded row is the same value at every pass and the driver folds it once.

## 3. Where `q(o)` comes from on each fold, and the one thing the design did not settle

Three sources, and the middle one is a decision this step took rather than inherited.

| when | what the row is scored against |
|---|---|
| no fraction fitted anywhere | no mixture at all — the plain formula of §3.3 |
| before the first pass | **flat**: every candidate allele equally likely in the contaminating population |
| every pass, and once more against the settled frequencies | the loop's own expected copies, summed over the sample's batch, with that sample's own copies taken out |

**⚑ The initialisation assembly's `q(o)` is not in the spec, and it is the one open choice here.**
§2's pseudocode builds the table once, before the initialisation pass, and §3 rules that that pass
scores the reads with *no prior at all*. With contamination on, the table it reads still needs a
`q(o)`, and the loop has produced no estimate yet — the per-sample expected copies are still the
scratch's `NaN` sentinel, which the copy fill refuses by name.

**It scores the reads alone: §3.3's formula, which is what this model computes wherever `c` is
zero.** That is the answer §3 asks the initialisation pass for — the reads, and nothing else — and
it needs no rule the model does not already have.

**⚑ The first version of this step used a flat `q(o)` instead, and the review's argument against it
is arithmetic on the model rather than taste.** A uniform distribution over candidate alleles is not
*saying nothing*: `c · q` is a floor under every observation's mixture that no genotype can lower, so
it compresses the differences between genotypes on the one pass whose whole purpose is to let them
speak. Computed on §3.6's own formula at a hom-ref genotype scoring four alternative reads, with
`ε̄ = 0.01`, a spread of 3 and `c = 0.05`: a flat `q = 0.5` makes those reads **28 Phred cheaper** to
explain than a converged `q = 0.05` does, against **3.7 Phred dearer** for scoring them with no
mixture at all. The reads-alone start is the closer of the two to where the loop settles.

**The third candidate was the prior's own seed**, and it goes out for §3's reason: the seed says a
locus is almost certainly invariant, which is the pull §3 rejects. **This is a decision this step
took, not a ruling** — it changes where every contaminated locus starts iterating and nowhere else.

**And an assembly that the plan's entry does not mention: one more after the loop stops.** The final
pass scores every sample against the settled frequencies, and with contamination on the table itself
depends on them — so it is assembled once more between the loop and the final pass. Without it the
final pass would score against the estimate the *last* pass started from. This is why the assembly
count is `passes + 2` rather than `passes + 1`.

**What it costs is worth naming**: the genotypes and the site quality then come from a table one
assembly *newer* than the one the convergence test looked at. At a locus that settled the two differ
by less than the convergence threshold — measured, the confidences move by less than a thousandth of
a Phred; at one that hit the pass cap the difference is whatever the last pass moved, and such a
locus is emitted flagged, which is what that flag is for.

## 4. `SequencingBatches`, and the rule nobody had settled

`arch/parameter_prepass_joint_fit.md` §1.6 specifies it and nothing had built it. It is now
`src/ng/parameter_estimation/joint/sequencing_batches.rs`, holding **two dense views of one
partition**: `BatchOfEachReadGroup`, which says which row of the frequency table a library's reads
are scored against, and `BatchOfEachSample`, which says which batch a sample's expected copies are
added into. The two are different lengths whenever a sample has more than one library, and identical
in length at one library per sample — which is every sample of every benchmark cohort here, and
exactly why they are two types.

**Deviation from the architecture, recorded rather than escalated:** the refusals are a
`SequencingBatchError` of this module's own rather than new variants on `JointFitError`. The type is
declared by the user and consumed by the caller; nothing in the fit produces it, and widening a
twenty-variant enum for it would have made the fit's error type answer a question the fit never
asks. The variant the architecture names, `ReadGroupNotBatched`, keeps its name.

**⚖ The rule for a sample whose libraries ran in different batches: refused, loudly.** The read
likelihood deliberately declined to invent one (C2's report §6), and there are three candidates —
pick a majority batch, average the two populations, or refuse. This refuses, naming the sample and
its batches. **Under the shipped default it cannot arise**, since one batch holds everything, so
what the refusal costs is a run that declares a batching splitting a sample and what it buys is that
nobody discovers the rule by reading genotypes. *This follows the owner's stated recommendation and
is not yet a ruling.*

## 5. The scatter, and the batch nobody was callable in

The plan asks for two buffers on `CallingScratch`, both `batches × alleles`. There are **three**, and
the third has a reason.

`fill_batch_allele_copies` takes the copies on one axis and the batching on the same axis, and
refuses a declared batch that no sample of that axis ran in — by name, as a batching that does not
describe the run. The loop's copies are on the **locus's** axis: one row per sample the locus is
called on, with the samples the candidate step ruled uncallable left out entirely. At a locus where
every sample of one batch is uncallable, summing over rows leaves that batch's row unwritten, and
the run dies on a check about the *batching* for a fact about the *locus*.

So the copies are scattered back onto the run's sample axis first — `expected_copies_by_run_sample`,
zero at every sample the locus has no row for. **Zero is the right answer and not a placeholder**:
an uncallable sample has no genotype estimate to contribute, which is exactly what the M-step
already does with it. `a_batch_with_no_callable_sample_here_is_a_row_of_zeros_rather_than_a_refusal`
is the fixture.

## 6. One performance change that the wiring forced, and one it did not

**The checks between the fractions and the batching are `read groups × batches`**, and the row loop
would have made them once a sample once a pass. Both sides of that, computed: at a thousand
libraries on four plates the checks are about 5,000 operations, so once a sample once a pass over a
thousand samples and seven passes is about **3.5 × 10⁷ at one locus** — against the row assembly's
own `samples × observations × genotypes × passes`, about **1.5 × 10⁶** at ten observations and a
six-allele diploid. Roughly twenty times the work the assembly does.

They are now made once per assembly, in a new `FrozenContamination`, which holds the run's half of
the mixture and hands out a `ContaminationMixture` over each refilled frequency table.
`ContaminationMixture::new` keeps its signature and is the two composed, so every test against it is
unchanged. `one_frozen_half_serves_two_frequency_tables` pins the reuse.

**What is still per row is named rather than fixed.** `with_frequencies` range-checks the table it
is handed, `batches × alleles`, and `fill_contaminant_allele_frequencies` rewrites the whole table
for every sample although only that sample's own batch row differs from the previous sample's. Under
the shipped default that is one batch — an allele's worth of work — and it stays small for a
plate-sized batching; it becomes the dominant cost only where a run declares roughly as many batches
as it has samples, which nothing produces today. **Closing it means changing what an earlier step's
function writes, not what this one checks**, so it is recorded rather than done.

## 7. Measurements

Every figure below was produced by the command beside it, in the container, on this tree.

| what | number | where |
|---|---|---|
| library target | **4,733 → 4,786** passing, 0 failed, 14 ignored | `cargo test --lib` |
| `ng::calling` alone | 707 → **738** passing | `cargo test --lib ng::calling` |
| release-held checks | **725** passing, 0 failed | `cargo test --release --lib ng::calling --all-features` |
| the loop's allocations, now on both paths | 1 passing | `cargo test --test ng_calling_loop_allocation --features dhat-heap` |
| broken intra-doc links | **28, all pre-existing and none on a line this change touched** | `cargo doc --no-deps --lib` |

**The passes a contaminated locus takes, measured on this step's own fixture** (three diploid
samples, three alleles, four reads of each allele at each sample, one library apiece, `c = 0.05`):
**7 passes against the same evidence's 4 with no fraction fitted.** `q(o)` moves between passes as
well as the frequencies do, so the loop has two moving quantities to settle rather than one. The
fold count is asserted as `passes + 2` at both a cap of 2 and uncapped, which is what makes
"independent of the pass count" a claim about the *emission* count rather than about a fixture.

**The counter's own reset is load-bearing, measured** (2026-08-26, `prepare_for_locus`'s
`emission_cost` reset deleted): the second locus called on one worker's scratch reports every field
doubled — `emission_builds: 2, emission_row_fills: 6, emission_evaluations: 36, table_assemblies: 2,
row_assemblies: 6` against `1, 3, 18, 1, 3`. A fresh scratch per locus hides it completely.

**What the last assembly is worth, measured at a locus that did not converge.** At a converged locus
the confidences move by less than a thousandth of a Phred whether or not it happens. Stopped after
one pass, dropping it moves the site quality from **119.720 to 119.742** and one sample's confidence
from **30.194 to 30.180** — small, and well above the `1e-4` the fixture allows.

**And what an allocation inside a contaminated pass costs, measured**: one `Vec::with_capacity`
beside the per-pass assembly takes the counted allocator from **8 blocks to 11**, one per extra pass;
reallocating the contaminant frequency table once a row once a pass takes it from **8 to 24** while
leaving every pointer fingerprint identical, because a freed block of the same size usually comes
back at the same address.

## 8. The release-held assertion battery

Every assertion this change adds outside a test module was downgraded to `debug_assert!` in one run
and the suite re-run under `--release`, where those checks vanish. **27 checks, and every one is
reached by a test that fails without it** — **29 tests fail**, across `ng::calling` and
`ng::parameter_estimation::joint::sequencing_batches`. (The `--release` lib suite has 8 failures of
its own on this tree, every one a `debug_assert!`-backed test in a module this change does not
touch — `genetics`, `ng::alignment`, `sample_summary`, `ssr`, `var_calling`. `ng::calling` itself is
0 failed in release, 725 passed.)

**Three checks are `debug_assert!` rather than release-held, each because no test can reach it.**
`FrozenParameters::gather` compares the contamination views against the batching in both directions
— a batching without views is a value nothing reads, views without a batching leave the mixture with
no row to score against — and `gather` is private with two callers, one always passing a batching and
one always passing none. `checked_axes`'s *read groups but no samples* is unreachable for a fact
about `ReadGroups`: it groups by sample, so a non-empty read-group table has at least one, and a test
says so rather than leaving the claim on trust. **A release check no test can reach is one the suite
cannot keep honest.**

**The two that a first battery missed are worth recording**, because the method is what caught them.
Downgrading *all* the checks at once hides a pair that shadow each other: `ContaminationMixture::new`
refuses a frequency table that is not a whole number of batches, and `with_frequencies` refuses one
that is not this run's batches by this locus's alleles, and with both gone the other's test still
failed on something. Isolating each and downgrading it alone is what shows which is reached — and
`with_frequencies`'s is the check the *loop* depends on, since the loop knows the batch count and the
one-step door derives it.

## 9. What the review changed, beyond the prose

**Seven deliberate defects survived all 4,776 tests when the reviews began**, and five were one
accident: every contamination fixture gave one library to one sample, which makes the two batchings
the same slice *and* makes a scratch row's index equal its sample's. The copies scattered onto the
row instead of the sample, the copies read from the sample instead of the row, the leave-one-out
subtraction reading sample 0's copies for everyone, and the two batchings swapped in either
direction — all passed. Three fixtures were rebuilt with the accident removed: an uncallable sample
placed **first**, four samples with genuinely distinct copies, and a run where one sample has two
libraries on one plate and another has one on a second. **And `FrozenParameters` now answers the two
batchings through `batch_of_sample(usize)` and `batch_of_read_group(ReadGroupId)`, so the
transposition is a type error rather than a number.**

**The three contaminant buffers were added to the allocation invariant and nothing exercised them** —
both readers of the pointer fingerprint and the counted-allocation test all called uncontaminated
runs, where the three are empty and every fingerprint is `(dangling, 0)` on both sides. Both halves
now have a contaminated arm. Two things came out of writing them: **two scratch buffers legitimately
exchange pointers with every pass** (the M-step swaps them rather than copying), so the fingerprint
list's order depends on the parity of the pass count and the new test compares the set; and **the
fingerprint cannot see a reallocation that returns to the same address**, which is what the counted
half is for.

**`row_assemblies` could not disagree with `table_assemblies × rows`**, because both were charged
from one argument outside the row loop. It is charged one row at a time now, so an assembly that
stops a row short moves it.

**The per-pass assembly's *effect* was pinned by no test** — only its cost. A golden-value fixture now
pins what a contaminated locus answers, at a converged locus and at a capped one, because at a
converged locus the last assembly moves the answer by less than the convergence threshold and only
the capped case can see it.

## 10. What this step did not do

- **The repeat-tract path is still refused at the seam's front door.** Its row needs a scoring
  context per `(read group, candidate)` whose outlier weight and reachable-length buffer no caller
  supplies; that is a step of its own, and the refusal's message no longer blames a parameter
  `FrozenParameters` now carries.
- **Nothing reports the fraction a run used.** That is E2b, and `SequencingBatches::is_default` —
  the one thing that can tell a declared batching from an assumed one, since the dense views cannot
  — is built and reachable from `RunParameters::sequencing_batches()` for it.
- **No CLI declares a batching.** `SequencingBatches::declared` exists and is tested; wiring
  `--sequenced-together` belongs with the run assembly.
