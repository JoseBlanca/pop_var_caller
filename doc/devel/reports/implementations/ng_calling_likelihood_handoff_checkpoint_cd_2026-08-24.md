# ng read likelihoods — handoff at Checkpoint C/D

*2026-08-24. Branch `ng-calling-likelihoods`, worktree `../pop_var_caller-calling-likelihoods`,
plan [`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). Written for
whoever picks this up next. It replaces
[the Milestone B handoff](ng_calling_likelihood_handoff_2026-08-24.md), which is now a historical
record.*

## Where the plan stands

**Seven of fifteen steps are done: A1, A2, B1, B2, C1, C2, D1.** That is **the whole generic
(SNP/indel) path** — plain, contaminated and censored evidence all scored. The next step is E1,
which begins the STR path and touches a different file.

| milestone | steps | state |
|---|---|---|
| A — module scaffold, shared vocabulary | A1, A2 | ✅ |
| B — the SNP/indel closed form | B1, B2 | ✅ |
| C — the contamination mixture | C1, C2 | ✅ |
| D — partial observations | D1 | ✅ **Checkpoint C/D reached** |
| E — the stutter distribution's three changes | E1, E2, E3 | ☐ next |
| F — the STR emission seam | F1, F2, F3 | ☐ |
| G — the censored term | G1 | ☐ |
| H — the STR row | H1, H2 | ☐ |

Suite: **4,354 passed, 0 failed, 14 ignored**; **162 in `ng::calling::likelihood`**, against 81
when this run began. Verified in the container and natively on the host.

## The eight things a successor must not rediscover

**1. An allele is the whole locus as a sample carrying it has it.** Not the reference with gaps
punched in. This is the single correction that shaped Milestone D, and getting it wrong sends you
looking for a map from an allele's bytes to the locus's positions — which does not exist, because
the merge computes each varying region's offset and discards it. It is not needed: a read from a
carrier shows the **start or the end** of that carrier's own sequence, so a partial is checked
against the allele's prefix or suffix. Spec §5.3 and arch §3 said otherwise until D1 corrected
them; if you find a document still promising a *gather*, it is stale.

**2. The witness counts locus positions; the bases are read content.** They differ by whatever
indel the read carried, so **the witness may never index the bases**. `PartialObservation`'s own
field docs say it. Every wrong version of D1 came from forgetting it.

**3. A partial constrains an allele only when every base it showed belongs to a run anchored at a
border.** Four shapes satisfy that (one run flush left, flush right, or flush both; two runs flush
both). Everything else leaves a run anchored to neither border whose bases absorb any
disagreement, so the test is **vacuous, not weak** — and D1 shipped a version that treated those
as contiguous, which *inverted* the verdicts: 14 nats charged against the allele the read actually
agreed with. `allele_is_compatible_with_partial` carries the argument.

**4. There is no ceiling on what a read is charged for being wrong.** Production pairs its
`1e-12` floor with a cap at a half; ng adopted both at A2 and the cap had to go, because a cap
binds on a single read and not on the fold of that read with others — a non-linear function of a
per-read quality, which spec §2.3 forbids outright. `MAX_BASE_ERROR` and the capped method are
deleted; `what_the_row_charges_a_poor_read_is_not_capped` is what stops them coming back.

**5. The error-spread table stores `m`, not `log m`.** B2 chose the logarithm with a sound
argument — §3.3's closed form takes no logarithm inside its loop — and C1 voided it, because
§3.6's mixture takes one there by specification and needs `m` itself. Arch §3 carries both
arguments; do not re-derive the first and reverse it back.

**6. `q(o)` leaves the sample out of its own batch.** The contaminant is somebody else by
definition, so the frequency it is drawn against excludes the individual being scored — the same
subtraction `fill_sample_concentration` already makes for the prior, sharing its
`COUNT_PATH_DESYNC_THRESHOLD`. `fill_batch_allele_copies` sums per locus;
`fill_contaminant_allele_frequencies` subtracts and normalises **per sample**. A sample alone in
its batch therefore gets the reference, which is right: no neighbours, no contaminating
population.

**7. `SequencingBatches` is specified and still unbuilt**
([`parameter_prepass_joint_fit.md`](../../ng/arch/parameter_prepass_joint_fit.md) §1.6). The loop
will need two views of it — `BatchOfEachReadGroup` and `BatchOfEachSample`, both in `ng::types` —
and **whoever builds it owes the rule for a sample whose libraries ran in different batches**. The
default (every read group together) is a complete answer, not a stub.
[`calling_loop.md`](../../ng/impl_plan/calling_loop.md) E2a and E2b own it, along with spec §3.6's
requirement that the run report the fraction it used per sample.

**8. `cargo clippy --all-targets --all-features -- -D warnings` is red on `main`**, in
`examples/ng_duplicated_class_harness.rs` and `benches/freebayes_bookkeeping.rs`. Validate with
`--lib --all-features --tests`, which is clean, and do not go chasing it.

## What Milestone E is, and why it is a change of gear

E1–E3 touch [`alignment/stutter.rs`](../../../../src/ng/alignment/stutter.rs) — **ng code, not
frozen production**, currently 15 tests — and the distribution there is **reused, not duplicated**:
one implementation with two consumers. E1 is a rename to the spec's vocabulary (`in_up`/`out_geom`
carry *in frame / out of frame*, which spec §1.3 bans) plus two doc repointings; it is mechanical
and its existing tests should stay green unchanged. E2 splits `MAX_SLIP` into two named cutoffs
and makes the truncated mass reported rather than silent. E3 adds the sums-to-one tripwire, which
spec §12's fourth test calls out as catching three silent failures at once.

**Everything from F onward is the STR path** and reads a different half of the spec (§4). The
generic path's types are done and should not need to change for it.

## Two habits this run paid for, repeatedly

**The review fan-out earns its cost, and what it found was never style.** Across seven steps it
found: a `NaN` that came back as the *most confident* read the model admits, because `f64::max`
returns the other operand; a spread table from the wrong ploidy truncating the genotype walk in
silence and leaving an unscored genotype the winner; a batch nobody was sequenced in coming back
indistinguishable from a batch that showed nothing; a divisor no fixture could see because every
sample in every fixture carried exactly two copies; and D1's inverted verdict. **Five of those are
wrong genotypes with nothing crashing.** Run it on every step, in worktrees, with mutation
testing.

**A fixture that agrees with your implementation is not the same as a fixture that pins the
property.** This bit twice in one milestone. The zero-contamination sweep claimed to catch a
reintroduced ceiling and could not — every profile's *fold* sat below the cap however poor its
single reads were. D1's fixtures could not tell the prefix/suffix rule from the positional gather,
because every one of them spelled alleles the reference's length and gave every read as many bases
as positions, which is exactly when the two agree. **Ask of every assertion which *wrong*
implementations also satisfy it**, and where a rule was chosen over an alternative, write the
fixture that separates them.

And the standing one: **every number about your own work is measured before it is written.** On
this run, wrong-first-time included a 5.51 that was 12.34, a 3 × 10⁻⁷ that was 3 × 10⁻⁵, a "seven
orders" that was nine, and a 1.69 that was 0.038. Every one was caught by asserting it in a test
rather than stating it in prose. Do that.
