# ng — a locus nobody can be called at is counted, not fatal

**Date:** 2026-09-01. **Branch:** `main`. **Owner's ruling**, taken at Checkpoint D: *"that
locus should be counted and reported at the end, but it shouldn't be an error that stops the
whole calling."* **Modules:** `src/ng/calling/mod.rs`,
`src/ng/calling/inference/summarise_condition.rs`, `src/ng/run/callers.rs`,
`examples/ng_call_cohort_end_to_end.rs`.

Not a step of any plan. It comes out of Checkpoint D, where D1's wiring made a panic reachable
from real data for the first time.

---

## What was wrong

The candidate step caps how many alleles a locus is called over, and a sample that had reads on
a sequence the cap cut is ruled **uncallable** there. `doc/devel/ng/spec/candidate_alleles.md`
§4.1 justifies cutting a sequence rather than refusing a locus on the ground that *most samples
stay callable* — and says nothing about the case where none does.

Until this change, that case hit an `assert!` inside `LocusGenotyper::call_locus`. So **one hard
locus ended a whole cohort's run**, with a Rust panic and a backtrace note.
`AlignedFilesVariantCaller::call_cohort` is what made it reachable from real data: everything
before D1 stopped at the evidence.

## What it does now

- **`LocusEvidence::callable_sample_count`** is the one place the question is asked. On a repeat
  tract it is the run's whole sample count, because a tract sets no sample aside.
- **`call_one_generic_locus` asks it before offering the locus to a genotyper**, and returns
  `None` where the answer is zero. The cohort observation is dropped and no record is made.
- **`CalledCohort::loci_with_nobody_to_call`** carries the ground of those loci, in genome
  order, so a run reports **where** and not merely how many. A non-empty list is worth acting
  on: raising `max_candidate_alleles` keeps more of what those loci vary over.
- **The assertion stays, reworded as a precondition on the trait** rather than as a claim about
  the data. It is now unreachable from every production caller — the review checked all 21 call
  sites and the other 20 are tests — and it stays for the SSR driver and the pool worker that
  will call `call_locus` next.

**It is a third fact, not a rename of either of the others.** *Too wide* is ground the merge
declined to assemble; *too quiet* is ground no sample varied at, which is counted nowhere by
design; this is ground that **was** assembled and where the cap then left nobody callable. And
it is not the per-sample `SampleGenotypeCall::Missing`, which is one sample set aside at a locus
other samples are still called at.

## ⚑ The trigger is narrower than it first reads, and the first draft's prose said otherwise

**`callable_sample_count` counts the run's samples, not the locus's covering ones**, and a
sample that covered nothing is **callable** — its evidence is empty, an empty sum is zero, every
genotype scores alike and the prior decides alone. So the guard fires only where **every sample
of the run** covered the locus *and* every one of them lost a sequence its own reads had earned.
In a 63-accession cohort at about three reads a position that is close to never, because some
sample almost always covers nothing.

**That is the right condition** — it is exactly the genotyper's precondition, whose scratch
cannot be prepared for no rows — but five doc comments described it as *every covering sample*,
which is a different and much commoner event. The review measured it: the two-alternative cohort
plus one sample whose reads lie on another contig produces **no** nobody-callable loci at all.
Corrected in all five places.

## Verification

| check | result |
|---|---|
| `cargo test --lib` | 5,818 passed, 13 ignored (5,813 before) |
| `cargo test --lib ng::run` | 378 passed (373 before) |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |

**Five tests, and the first of them reaches the old panic.** Removing the new guard and
re-running it panics at `summarise_condition.rs` with *"every one of the 2 samples at contig
0:15-15 was ruled uncallable"* — so the fixture provokes the case rather than describing it.

The fixture is two samples that **each** show `C` on two reads, `G` on two and the reference on
two, at a cap of one alternative: both alternatives clear the merge's floor, both are earned by
both samples, and whichever the cap cuts, every sample of the run has lost a sequence its own
reads earned. The two tests beside it are what say the first is about the cap rather than about
the reads — the same cohort at the shipped cap of six calls that locus over three alleles, and a
cohort where only one of three samples loses its allele is still called with that one sample set
aside.

**On real data at the shipped cap it does not fire.** Six tomato accessions over 400 kb of SL4.0:
8,411 loci called, **0 where the cap left nobody callable**. That is a measurement rather than a
reassurance — the number is printed by the end-to-end probe on every run.

## What the review changed

A 14-mutation pass left five survivors. Two were the fixture rather than a defect and three are
now dead:

- **A guard written `<= 1` rather than `== 0` survived everything**, and it would empty a
  **single-sample run**: at one sample, "fewer than two callable" is every locus. That is the
  thinnest end of the range this caller commits to and no other fixture in the file exercises
  it. There is now a one-sample test.
- **Nothing asserted the three lists are disjoint**, which is what the type's own documentation
  claims: a locus counted in two would be reported twice and make the assembled total wrong in
  both directions.
- **Nothing pinned "in genome order"**, because one locus cannot tell an ordered list from a
  reversed one. The fixture now places both alternatives at two positions, so a cohort of it
  has two such loci.

The two that remain are recorded rather than forced: a one-base span cannot show that the
recorded region's end is right — the shared fixture reference is a hundred `A`s, so no
observation spans more than one base, which is the same limitation already recorded against
`loci_too_wide_to_assemble`; and the probe's own output line is compiled by `clippy
--all-targets` but exercised by no test, since `cargo test --lib` never builds examples.

**And six prose claims were wrong**, beyond the subject error above: the wrong document was
cited for the truncate-never-refuse ruling (it is `spec/candidate_alleles.md` §4.1, not the
architecture's), a paragraph about the allele table was duplicated, `CalledCohort`'s summary
still said "three" of four fields, an `expect` message claimed the oracle disagreed about
callability when the oracle is never asked, and **the probe told an operator to raise
`--max-candidate-alleles`, a flag that does not exist** — the command surface is Milestone F.
The probe now also prints the loci the merge assembled, because `called_loci` alone stopped
being that number.

## ⚑ A correction owed to the spec, which is the owner's to make

`doc/devel/ng/spec/candidate_alleles.md` §4.1 rules that selection cuts an allele rather than
refusing a locus, *because most samples stay callable*. It does not cover the case where none
does, and the code now has an answer for it: count the locus, report its ground, call nobody
there. Recorded here and in `PROJECT_STATUS.md` rather than written into the document.
