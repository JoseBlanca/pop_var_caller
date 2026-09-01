# ng direct mode, step D2 — the sample-order join

**Date:** 2026-09-01. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step D2. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md) §5.1.
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §3.2 —
*"this is where the three sample numberings meet"*. **Module:** `src/ng/run/callers.rs`.

**Tests and fixtures only.** No production code changed: the join was built with the calling
loop and with D1's wiring, and what it lacked was a cohort that could tell a correct join from a
wrong one. Its own commit because the plan asks for one, and the reason is that a defect here
produces wrong genotypes rather than a crash.

---

## The three numberings, and what a fixture has to do to separate them

- **The merge's** holds only the samples that covered the locus, each carrying its own index in
  the run's order.
- **The run's** is `ReadGroups::read_groups_per_sample`'s first-seen order. Every per-sample
  list the calling loop is *given* is in it — the evidence, and the model parameters.
- **The calling scratch's rows** are the run's samples with the uncallable ones closed up, and
  the loop's own working buffers are indexed by these.

**A permutation of any of the three produces a well-formed answer for every sample**, which is
why nothing about a genotype reveals one. What a swap changes is *which* answer, and only where
the samples are scored under something of their own. Here that is each sample's inbreeding
coefficient — how much the prior expects homozygotes — declared **by name** and resolved
against the run's order.

## The cohort

| run index | sample | at `chr1:15` | merge entry | scratch row |
|---|---|---|---|---|
| 0 | `zeta` | two reads show `C` | 0 | 0 |
| 1 | `nu` | two reads show `G` | 1 | — set aside |
| 2 | `alpha` | no reads in the analysed ground | — | 1 |
| 3 | `mu` | two reads show `C` | 2 | 2 |

**`alpha` separates the merge's numbering from the run's** by covering nothing where the others
cover. **`nu` separates the scratch's rows from both**: its `G` is the cohort's lower-ranked
alternative — two reads against the four behind `C` — so at a candidate cap of one alternative,
selection cuts it and rules `nu` uncallable for having earned a sequence the cap removed. `mu`
is then the run's sample 3, the merge's entry 2 and the scratch's row 2.

**`zeta` and `mu` bring identical reads**, so any difference between their calls is a difference
in what they were scored under and nothing else. **The names are not alphabetical**, and the
hazard that guards against is real rather than theoretical: `DeclaredInbreeding` holds its
per-sample values in a `BTreeMap`, so a defect that zipped that map's key order onto the run's
samples would hand `zeta` what `mu` was declared.

**That the fixture is three different numberings is asserted, not described.** Without that
assertion every other test here could be passing on the identity permutation, where a run that
indexed one list by another would look correct.

## What the discriminator is, and what it is not

**At this depth the coefficient moves the confidence, not the genotype.** At the fixture's base
quality of 30, two alternative reads cost a homozygous-reference genotype about a millionth, so
the reads decide the heterozygote and the prior only moves how sure the caller is of it.
Measured on this cohort:

| sample | declared | call |
|---|---|---|
| `zeta` | outbred | `0/1` at **55.449852** Phred |
| `mu` | nearly fully inbred | `0/1` at **33.363934** Phred |
| `alpha` | nothing, and no reads of its own | `0/1` at **2.2188103** Phred |
| `nu` | nothing; the cap cut its allele | no genotype |

So the tests compare the **whole call**. A genotype comparison alone would have been blind to a
parameters list joined by the merge's entry — measured, that defect leaves both samples at `0/1`
and 55.450. It would **not** have been blind to a wrongly joined *evidence* list, which moves
the genotypes, and the genotype assertions are what guard that.

**The oracle** is `swapping_two_samples_coefficients_swaps_their_calls`: the same four files,
reads, order and cap, with only which sample each coefficient names exchanged, requiring the two
calls to exchange exactly — and the two samples between them, named by neither declaration, to
not move at all.

## Verification

| check | result |
|---|---|
| `cargo test --lib` | 5,813 passed, 13 ignored (5,809 before this step) |
| `cargo test --lib ng::run` | 373 passed (369 before) |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |

**Four mutations at the joins, and D2's own four tests catch every one.** Each was injected,
run, and the tree restored and re-checked:

| mutation | D2's tests alone |
|---|---|
| index the run's per-sample rows by the merge's entry | all four fail |
| read selection's leftovers by the run's sample rather than the merge's entry | all four fail |
| shrink the run's sample count | all four fail |
| sort the run's sample names | one fails |
| read a sample's coefficient by its scratch row, not its run sample | two fail |
| emit the per-sample evidence views in reverse run order | two fail |

**Scoping the run matters and is the reason these were run twice.** The first two mutations are
inside `evidence_shaping`, which has a substantial suite of its own, so `cargo test --lib ng::`
would have caught them whatever D2 did. Run against `ng::run::callers::the_sample_order_join`
alone, they still die — which is the claim this step is entitled to make.

**And the two the review found in the deeper code die here too**, which the first fixture could
not have managed. The correctness review wrote fourteen mutations against the first draft and
two survived everything in the crate:

- **the calling loop reading each sample's coefficient by its scratch row rather than by its
  run sample** (`summarise_condition.rs`). On a cohort where every sample is callable the two
  indices are equal, so nothing anywhere in `ng::` could see it. On this cohort `nu` has no row,
  so `mu` is run sample 3 and row 2 — and the mutation now fails two of D2's four tests.
- **the per-sample evidence views emitted in reverse run order** (`evidence_shaping::fill_views`).
  It survived because `zeta` and `mu` were given identical evidence on purpose, and with three
  samples the middle one is a fixed point under reversal. With four, reversing exchanges `nu`'s
  reads with `alpha`'s — and `nu` is the sample the cap sets aside, so the pattern of which
  samples are called moves. Two of the four tests fail.

## What the review changed

Two reviews ran in parallel over the step's diff, each in its own worktree.

### The fixture met two thirds of the oracle

**The first draft separated the merge's numbering from the run's and left the scratch's rows as
the identity.** A sample becomes uncallable only when the allele *cap* cuts a sequence its own
reads earned; the first fixture had one alternative against a cap of six, so nothing was cut,
every sample was callable, and a row index equalled a run index. Of the three pairs the plan's
oracle asks about, exactly one was separated — and nothing said so, so a later reader had no way
to see it. `nu` and the narrowed cap are the fix.

### A claim of this step's own was wrong by about five orders of magnitude

The fixture's documentation said two reads each way "say a heterozygote and a homozygote almost
equally loudly, so what decides the call is the prior". Measured, the caller gives the
heterozygote about 350,000:1 at the outbred end and still 2,200:1 at the inbred end: the
likelihood decides the genotype and the prior moves only the confidence — which is what the
module header said twelve hundred lines further down. The same wrong sentence appeared in the
fixture door's own doc.

### The walk tallies were checked at one index

`a_sample_that_covered_nothing_is_called_by_the_prior_alone` asserted a single name. Index 2 of
four is a fixed point under reversal about the middle, and under a zip that drops one sample and
leaves the rest shifted the assertion holds on a list wrong in both length and pairing. It now
asserts the whole vector.

### And four smaller things

- **An ordering assertion written on `Option<f32>`.** `None < Some(_)`, so a regression that
  made the inbred sample `Missing` would have passed the assertion while its message printed
  `None` where it claimed a Phred score.
- **A duplicated genotyper helper**, doc comment and all, whose copy pointed at the original —
  two things to keep in step the first time the shipped emission model or prior changes.
- **A fixture door that bypassed the one function every other door funnels through**,
  re-spelling the whole constructor and hardcoding five settings. A new field on
  `AlignmentInputs` is a compile error at every door; a change to *which* settings the fixtures
  open with is not — and that is exactly the drift that produced D1's oracle regression.
- **Two clauses that announced importance rather than adding a fact**, and a `genotype_of` that
  returned a `String` where `Genotype` already compares and prints.

### One thing the reviews raised that is the owner's, not this step's

**The plan's D2 oracle says "swapping any two changes a called *genotype*".** On a fixture at
real read depth it changes the *call* and not the genotype, because the reads decide the
genotype long before the prior does. Making the genotype itself flip would need evidence
contrived to be ambiguous — low base qualities at the variant position — which would be a
fixture built to satisfy a sentence rather than to resemble a run. The code has overtaken the
plan's wording; the wording is the owner's.

## What this step does not do

- **It does not pin the join for repeat tracts.** A tract sets no sample aside, so its scratch
  rows are the run's samples exactly, and both tract generator slots are unfilled anyway.
- **It does not pin the read-group numbering**, which is the fourth numbering in the run and is
  keyed by read group rather than by sample. Nothing in this step moves it.
