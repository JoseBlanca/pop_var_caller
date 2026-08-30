# ng calling loop — E2a review: the contaminant frequency, per locus and per sample

**Step:** E2a of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md).
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.
**How it was run:** three agents, each in its own git worktree detached at `743aace3` with the
working tree's patch applied — one on the arithmetic and the control flow, one on the tests and the
cost invariants, one on design conformance and the accuracy of every claim.
**Verdict:** Request-changes. **1 Blocker / 11 Majors / 16 Minors / 9 Nits**, and **16 of 96 of the
diff's own claims wrong** — every one of them a mechanism rather than a number.

---

## 1. The headline, in one line each

**No arithmetic defect was found in the shipped code.** The split of the genotype-likelihood
table's build was compared term by term against the pre-change file and moves nothing; the driver's
ordering is right; no unwritten sentinel reaches a fill; no check was lost when
`ContaminationMixture::new` was split in two.

**What was wrong was the evidence.** Seven deliberate defects survived all 4,776 tests, and five of
them were the same accident wearing three faces.

**And the one Blocker was a sentence that had been true and was not any more**, standing in eleven
places — the fourth time this plan has shipped that shape.

---

## 2. The Blocker: the ruling reached the code and one document, and eleven statements still said
the opposite

The owner's ruling splits the build: the emission stays computed once per locus, the per-genotype
row is assembled again at every pass wherever contamination is on. The code does that. Three
surviving statements said the table is never rebuilt, and a grep for the rest of the sentence found
**eleven**, of which the worst was
[`summarise_condition.rs`](../../../../src/ng/calling/inference/summarise_condition.rs)'s own
`run_frequency_loop`: *"It reads the genotype likelihood table the scratch already holds and never
rebuilds it"* — five lines above the new parameter that rebuilds it.

Two were in the spec itself: §5's table of what moves between passes listed each sample's read
likelihoods as *no*, and §8's cost paragraph counted one build. **Both are now amended in place with
the ruling and its date**, because the ruling was already taken and what was missing was the record
of it, not a new decision.

## 3. The Majors that were tests unable to fail, and the accident behind them

**Every contamination fixture gave one library to one sample.** That makes
`BatchOfEachReadGroup` and `BatchOfEachSample` the same length *and* the same values, and it makes a
scratch row's index equal its sample's index. Three transpositions and two mis-indexings passed the
whole suite:

| the defect | why no fixture saw it |
|---|---|
| the copies scattered onto the **row** index instead of the sample the row names | the one fixture with an uncallable sample put it **last**, where the two indices coincide |
| the copies read **from** the run-sample index instead of the row | the same |
| the leave-one-out subtraction reading **sample 0's** copies for every sample | three of the fixture's four samples carried identical copies, and the asserted one matched sample 0 |
| the sample batching handed to the **read-group** half of the mixture | the two views are identical at one library per sample |
| the read-group batching handed to the **sample** half | the same |

**All five now die.** Three fixtures were rebuilt with the accident removed — an uncallable sample
placed *first*, four samples with genuinely distinct copies, and a run where one sample has two
libraries on one plate and another has one on a second, so the two views differ in length and in
content. **And one fix is structural rather than a test:** `FrozenParameters` now answers the two
batchings through `batch_of_sample(usize)` and `batch_of_read_group(ReadGroupId)`, so the
transposition is a type error. The fixture stays as the guard if the two are ever collapsed back
onto one call.

## 4. The three contaminant buffers were added to the allocation invariant, and nothing exercised
them

`buffer_fingerprints` gained the three buffers this step adds, with a doc saying they are *"exactly
the property this measures"*. Both of its readers, and the counted-allocation test beside it, called
uncontaminated runs — where all three are empty and every fingerprint is `(dangling, 0)` on both
sides. **This is the third time this plan has shipped an assertion that looks at a view which comes
back empty either way.**

Both halves now have a contaminated arm. Writing them turned up two things worth keeping:

- **Two scratch buffers legitimately exchange pointers with every pass** — the M-step swaps the
  cohort's expected copies with the previous pass's rather than copying either — so the fingerprint
  list's *order* depends on the parity of the pass count. The existing test passes only because both
  its runs take an even number of passes. Documented, and the new test compares the set.
- **The fingerprint cannot see a reallocation that returns to the same address**, which is what a
  freed block of the same size usually gets. Measured: reallocating the contaminant frequency table
  once a row once a pass leaves every fingerprint identical and takes the counted allocator from
  **8 blocks to 24**. That is what the counted half is for, and its doc now says so.

## 5. The design Major: the first assembly's `q(o)`, and the code changed

The step shipped a **flat** first guess — every candidate allele equally likely in the contaminating
population — for the assembly the prior-free initialisation pass reads, on the argument that it is
what §3 gives the genotype prior for the same reason.

**The review's counter-argument is right and is arithmetic on the model.** A uniform `q` is not the
neutral answer: `c · q` is a floor under every observation's mixture that no genotype can lower, so
it compresses the differences between genotypes on the one pass whose whole purpose is to let those
differences speak. Computed on §3.6's formula at a hom-ref genotype scoring four alternative reads,
`ε̄ = 0.01`, spread 3, `c = 0.05`: a flat `q = 0.5` makes those reads **28 Phred cheaper** to explain
than a converged `q = 0.05` does, against 3.7 Phred *dearer* for scoring them with no mixture at all.

**The initialisation assembly now scores the reads alone** — §3.3's formula, which is what this model
computes wherever `c` is zero. It is the closer of the two to where the loop settles, it needs no
rule the model does not already have, and it deletes a code path. The test asserts the strong form:
that first table is *entry for entry* the table an uncontaminated run gets.

## 6. Two release-held checks no test could reach, and what happened to each

The battery downgrades every added assertion to `debug_assert!` and re-runs under `--release`. The
step declared two exceptions — the pairing checks inside `FrozenParameters::gather`, unreachable
because its only two callers each always pass one shape. **That claim held.** Two more were not
reached and were not declared:

- `FrozenContamination::with_frequencies`'s length check. Unreachable through
  `ContaminationMixture::new`, which derives the batch count *from* the table's length — but the
  loop uses the two-step door, which knows the batch count already and so catches the case the
  one-step door cannot: a table of the right shape for a different run. **It has a test now**, and
  downgrading it alone fails exactly that test.
- `checked_axes`'s "read groups but no samples". Genuinely unreachable: `ReadGroups` groups by
  sample, so a non-empty read-group table has at least one. **Demoted to `debug_assert!`**, with a
  test that says why it is unreachable rather than merely untested.

**Final battery: 27 added release-held checks, 29 tests fail when they are downgraded**, against a
baseline of 8 release failures in modules this change does not touch.

## 7. What the counters could not tell apart, and what now can

Three defects were caught by the cost counter and by nothing else: never assembling inside the loop,
skipping the assembly against the settled frequencies, and using the first guess at every pass. So
the *cost* of the per-pass assembly was pinned and its *effect* was not — the contaminated fixtures
asserted a direction and two counters, and a direction survives all three.

**A golden-value fixture now pins what a contaminated locus answers.** Measured: moving the
head-of-pass assembly to after the E-step's row loop leaves all 733 other `ng::calling` tests green
and moves the cohort's expected copies by about **3.3 parts in a million** with the pass count
unchanged. And because at a converged locus the last assembly moves the answer by less than the
convergence threshold — which is the bound the code documents — the fixture pins the **capped** case
too: stopped after one pass, dropping that assembly moves the site quality from **119.720 to
119.742** and one sample's confidence from **30.194 to 30.180**.

**And `row_assemblies` was a field that could not disagree with `table_assemblies × rows`**, because
both were charged from the same argument outside the row loop. It is now charged one row at a time,
so an assembly that stops a row short moves it.

## 8. The claim that was wrong about this project's own design

The step called the emission *"the expensive half"*, following the spec's word for it. Counted on
the module's own contaminated fixture — 3 samples × 3 observations × 3 alleles, diploid, so 6
genotypes and 9 assemblies — the assemblies cost about **486 multiply-adds and 54 logarithms**
against the emission fill's **9 charged-error calls**. The emission is genuinely the expensive half
for the byte comparison per `(partial read, candidate)`, and on the repeat-tract path for an
alignment per `(observation, candidate)`. **On plain SNP/indel evidence the assembly is the larger
side, and what the split buys there is the count rather than the clock** — that emission evaluations
cannot grow with the pass count. Both doc comments now say so.

## 9. Claims

**96 checked, 16 wrong**, split 35/61 across the two agents that counted. Every wrong one was a
mechanism; **every counted figure that was re-derived was right** — the reset measurement, the
`passes + 2` assembly count, `Σ_s` at 27 evaluations, the batch copy tables, the frequency rows, the
twenty per-locus buffers, and the tomato/SRA sentence quoted from the architecture.

The wrong ones, in one line each: the table is never rebuilt (eleven places); the assembly runs
"once per pass" where it runs `passes + 2` times; the flat first guess "says nothing"; the empty
contaminant tables "would not fail on their own" (the copy fill refuses them first, as a shape
mismatch — what the guard adds is naming which call was skipped); `is_default` "nothing has built"
(this step built it); "the same eight things" in a nine-argument constructor; a test doc naming the
wrong mechanism for what its fixture catches; the zero fill that was dead where it stood; "the
expensive checks are made once" where they are made once *a pass*; "the loop is where it belongs"
for a return value the loop then dropped; and the expensive-half framing of §8 above.

## 10. What was left undone, deliberately

- **The per-row cost of the contaminant block.** `with_frequencies` range-checks `batches × alleles`
  per row and the fill rewrites the whole table per row, though only one row of it differs between
  samples. Under the shipped default that is one batch — an allele's worth of work — and it stays
  small for a plate-sized batching; it dominates only where a run declares roughly as many batches
  as it has samples, which nothing produces today. Closing it means changing what
  `fill_contaminant_allele_frequencies` writes, which is an earlier step's function and an earlier
  step's contract. **Documented with both sides of the arithmetic rather than fixed.**
- **`SequencingBatches::declared` has no producer.** No command-line flag carries a batching, so
  every run today gets one batch holding everything. The module says so now, so the refusals are not
  mistaken for coverage.
- **The plan's E2a entry says two buffers on `CallingScratch` and there are three.** The third is
  argued in the implementation report and tested; the entry was not amended.
