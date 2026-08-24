# ng read likelihoods — C2: the contaminant is drawn from the neighbours, not from the cohort

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step C2 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone C, on
top of `d9028ff8`'s parent.*

## 1. What it is

C1 took `q(o)` as a parameter. This step produces it.

`q(o)` is **how often the samples sequenced beside this library carry the allele this read
showed**. A batch's frequency is its samples' expected allele copies added up and divided by
their total — the loop's own numbers, so nothing is estimated here that the loop has not
estimated already. That is what makes it a lookup rather than a fit, and what makes it right to
recompute every iteration instead of freezing it beside the contamination fraction (spec §3.6,
corrected 2026-08-24).

```text
                    Σ  expected copies of allele a in sample s
        q(a | b)  =  s ∈ b
                    ─────────────────────────────────────────────
                     Σ    Σ  expected copies of allele a′ in s
                    s ∈ b  a′
```

## 2. The shape, and the thing that made it fit

**The row needs no sample identity for any of this**, which is the reason the design stayed
small. Sequencing batches are declared over **read groups**, not samples — a sample's libraries
may have run on different flowcells, so the read group is the grain the header gives
([`parameter_prepass_joint_fit.md`](../../ng/arch/parameter_prepass_joint_fit.md) §1.6) — and
every observation the row walks already carries its read group. So:

- `ContaminationMixture` holds the frequencies **batch-major**, one row of `allele_count` per
  batch, beside a `BatchOfEachReadGroup` saying which batch each read group ran in.
- `contaminant_frequency_of(read_group, allele)` reads the row the batching points at.
- `genotype_log_likelihood_row`'s signature is unchanged, and `GenericSampleEvidence` still
  carries no sample field.

**The default batching is one batch holding the run** — `BatchId::ALL_TOGETHER` for every read
group —
so the table is a single row, every read group reads the cohort frequency, and nothing branches
on the batching's absence. A run that declares no batching loses nothing it had.

## 3. What the batch axis is worth, measured

The same three alternative reads at a summed log error of −21 and a 4% contamination fraction,
scored under a reference homozygote: once from a library that ran beside samples carrying the
alternative at 1 in 2, once from one that ran beside samples carrying it at 1 in 1,000.

**The two libraries differ by 12.34 nats — 54 Phred — at the same locus, on the same reads, at
the same fraction.** Per read that is `ln(0.0203 / 0.000332)`: in the first batch the
contaminant route `0.04 × 0.5` swamps the misread's `0.96 · e⁻⁷/3`, and in the second
`0.04 × 0.001` does not. **A row that ignored the batching would give both the same number**,
which is the whole of what this step added.

And the default is pinned rather than argued: under one batch the row is bit-identical whichever
read group an observation came from.

## 4. The floor, and why its size is an argument rather than a default

An allele a batch never showed gets `MIN_CONTAMINANT_FREQUENCY`, `1e-12`, not zero.

**The temptation is to make it statistical, and that would be the wrong direction to be wrong
in.** A batch that never showed an allele is not proof its samples lack it: with 63 diploids —
126 chromosomes — an allele seen zero times could still sit near 1 in 42. But encoding *that*
uncertainty here says "this read might well be contamination" at every candidate the cohort is
thin on, which is most of them, and spec §3.6 names exactly that failure as the one to watch:
contamination attributed too readily suppresses real heterozygotes by explaining their
alternative reads as somebody else's.

So the floor is **defensive**, and it is set where it cannot compete with the route beside it:
at a 3% fraction a floored allele contributes `3 × 10⁻¹⁴` against a misread's `3 × 10⁻⁵` at a
middling Phred 40 — **nine orders of magnitude below**, pinned by a test.

**Whether it should instead be a pseudocount over the batch's copies is a modelling decision
this step did not take.** It belongs in the producer, it has a genotype effect, and it is the
owner's.

## 5. A batch that shows nothing

A batch whose samples hold no copies at all — no coverage at this locus — has no frequency to
read off, and its zero total would give a row of `NaN`. A `NaN` reaches a logarithm and then an
argmax, where **every comparison against it is false**, so a genotype would be picked in silence.
It gets the reference allele and the floor elsewhere instead: the honest statement of what a
batch with no evidence has to say about a contaminant, and it keeps every row a distribution.

## 6. What is deliberately not here

**`SequencingBatches` is still unbuilt**, and building it belongs to the joint fit's plan rather
than this one. The mixture takes the batching as a `BatchOfEachReadGroup`, which that type produces
trivially when it lands.

**Which batch a sample belongs to when its libraries ran in different ones is a rule this step
does not invent.** Batches are over read groups; a frequency is over samples; so something has
to say what a sample split across two flowcells contributes to each. That rule belongs with
`SequencingBatches`, so the producer takes the sample's batch as an argument rather than deriving
it. Every sample of every benchmark cohort here has one library, so the case is not exercised by
anything in the tree — which is a reason to name it, not to assume it away.

**And the fraction used still has no route to the run's output.** Spec §3.6 requires it per
sample, because a genotype computed at `c = 0.03` and one at `c = 0` are otherwise
indistinguishable; no step of this plan owns it. Carried forward from C1's report unchanged.

## 7. What the reviews changed

Review agents ran on the committed step, each in its own worktree, over reliability, errors,
defaults, naming, idiomatic and smells.

**Two ways for a batching to be wrong and say nothing, both now refused.** The producer never
writes a batch no sample is in, so `out.fill(0.0)` left it at a zero total and it came out
through the no-evidence fallback — **a row indistinguishable from a batch that really was
sequenced and really showed nothing**; measured, a table sized for three batches against a
batching naming only the first returned `[0.999999999999, 1e-12]` for the other two. And the
mixture checked that every batch a read group *names* has a row, but not that every row is
named, so a three-row table with an all-default batching constructed cleanly, every read group
read row 0, and **a run that declared batches was scored against the cohort frequency**. The two
are mirrors, and without both a defaulted batching slipped past each end. The function's own
`# Panics` section had promised the first of these all along.

**A row of finite copies whose sum overflows.** Each slot is checked as it accumulates and the
total was not, so one sample holding `1e308` copies of two alleles gave every ratio as
`finite / inf` — zero — and every allele was lifted to the floor. The row came back finite,
plausible, and saying the neighbours carry nothing.

**Two batchings that a call site could not tell apart.** `fill_contaminant_allele_frequencies`
takes one keyed by sample and `ContaminationMixture::new` one keyed by read group; both were
`&[BatchId]`, so each function accepted the other's argument. **The length checks catch the swap
only when the sample count differs from the read-group count — and at one library per sample
they are equal**, which is every sample of every benchmark cohort here. Sample order and
read-group order are minted by different rules, so a silent mis-key is a wrong `q(o)` for every
observation, worth up to 12.34 nats a read by §3's own measurement. They are now
`BatchOfEachReadGroup` and `BatchOfEachSample`, on the same argument the allele-copy views in
`genotype_prior` already won: two types for one shape, because the swap is otherwise invisible.

**And `BatchId::ONLY` was a fifth word** for what the architecture already calls `all_together`,
beside the constant, the helper and two sentences of prose. It is `BatchId::ALL_TOGETHER` now.

Smaller: the fill returns **how many batches took the no-evidence fallback**, because otherwise
it is invisible — at one candidate allele a batch full of evidence and a batch with none are
bit-identical; the reference allele is spelled by its constant rather than a bare `0`; the row
loops read `chunks_exact`; the whole-table postcondition is `debug_assert!`, since it sweeps once
per locus to prove something no input can break; and the batching fixture, which had been copied
into both test modules, lives in one place.

**Recorded as checked, so it is not re-checked:** `allele_count = 1`, denormal totals, and a true
ratio of `1e-600` all behave; `BatchId` deliberately does not derive `Default`; the batching is a
required positional parameter rather than an option.

### What mutation testing added

**Sixteen mutations run, six survived, none changed no behaviour** — and the one that mattered
is a divisor. **No fixture separated *divide by the batch's own copies* from *give every sample
one vote***, because every sample in every fixture carried exactly two copies: `FOUR_DIPLOIDS` is
four rows summing to 2, and even the deliberately fractional one is `[1.5, 0.5]` and
`[0.25, 1.75]`. A divisor of two copies per sample passed all 127 tests. There is now a
mixed-ploidy fixture — a tetraploid at 3:1 beside two diploids, five reference copies and three
alternative out of eight — where the wrong divisor gives `[0.833, 0.5]`, **a row summing to
1.33**, and a property test over 3 to 20 samples, 1 to 5 alleles and 1 to 4 batches pinning the
two laws this producer obeys: every row is a distribution, and multiplying every copy by the same
number changes nothing. Both catch it, checked by putting the mutation back.

Three release-mode guards had no test and now do: a negative copy count (the finiteness check
alone does not catch it, and every other fixture reaching that guard uses a `NaN`), a frequency
buffer that is not a whole number of rows, and `BatchId`'s own rendering.

**And a batch of one sample is now written down.** `q(o)` comes back as that sample's own
genotype — a lone heterozygote gets `[0.5, 0.5]` — which is §8's open question, and a test saying
so is better than a gap. It pins the behaviour without endorsing it.

## 8. Two questions this step raises and does not answer

**A batch of one sample makes `q(o)` the sample's own genotype.** The frequency is summed over
the samples in the batch — **including the sample being scored** — so an alternative homozygote
alone in a batch has `c · q(alt)` explaining its own alternative reads. At 63 samples its own
contribution is 1 in 63 and the circularity is negligible; at one it is everything.

**The prior already solved this problem and `q(o)` did not inherit the solution.**
`fill_sample_concentration` takes the cohort's copies and subtracts the sample's own —
leave-one-out — precisely because a sample must not be its own evidence. The same argument
applies here and more directly: the contaminant *is* somebody else, by definition, so the sample
itself does not belong in the population its contaminant is drawn from. **The reason C2 did not
do it is that it changes the shape**: a leave-one-out `q` is per `(sample, batch, allele)`, and
the row deliberately has no sample identity — which is the whole reason this step stayed small.

Production takes the blunt route at this size: below five samples in a batch it zeroes the
contamination fraction rather than trusting the frequency. That is a third option and the
cheapest.

**Recommendation: subtract the sample's own copies, and pay for it in the loop rather than in the
row.** The loop knows which sample it is scoring, so it can fill a per-sample row of the
frequency table before calling the row — the table is `batches × alleles` and would become
`1 × alleles` per sample, refilled per sample. That keeps the row's signature and its
sample-blindness intact. But it is a modelling decision with a genotype effect at small batches,
so it is the owner's, and it should be settled before the loop plan fixes the calling order.

**The second question is the floor's size**, carried from §4: defensive at `1e-12`, or a
pseudocount over the batch's copies. The two are the same question at different ends — what a
batch that shows little or nothing is allowed to say about a contaminant.

## 9. Validation

All in the container, on the committed tree:

- `cargo fmt --all -- --check`: clean.
- `cargo clippy --lib --all-features --tests -- -D warnings`: clean. The repo-wide
  `--all-targets --all-features` run is red on `main`, in `examples/ng_duplicated_class_harness.rs`
  and `benches/freebayes_bookkeeping.rs` — pre-existing and out of scope.
- `cargo test`: **4,330 passed, 0 failed, 14 ignored**; 138 of them in
  `ng::calling::likelihood`, against 113 at C1.
