# ng read likelihoods — D1: a read that ran out is a prefix or a suffix of whatever produced it

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step D1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone D, on
top of `9d797b3e`'s parent. **This completes the generic path — Checkpoint C/D.**

## 1. What it is

Some reads run out inside a locus: they cover part of it and nothing outside. Until now the row
threw them away. D1 scores them.

## 2. The correction the step turned on

The plan, the architecture and spec §5.3 all described the rule as a **positional restriction** —
take the allele's projection, restrict it to the positions the read witnessed, compare. That
reading needs a map from an allele's bytes to the locus's positions, and for an allele that
inserts or deletes **no such map exists anywhere in the code**: the merge computes the offset of
each varying region and discards it, storing the allele as a flat byte string.

**The rule is not a positional restriction** (owner, 2026-08-24, correcting the question this step
asked). **An allele is the whole locus as a sample carrying it has it**, not the reference with
gaps punched in. So a read from a carrier shows the start or the end of that carrier's own
sequence, and the comparison is against the allele's **prefix** or **suffix** — nothing is
assembled and no map is needed. `WitnessedLocusPositions::is_flush_left` and `is_flush_right` are
documented as precisely those two constraints; this is their consumer.

**Three consequences.** There is no gather, so the buffer §5.3 said this step owed a home does not
exist — the step's scratch is the compatibility cache alone. Spec §5.3 and arch §3 are corrected,
because both stated the consequence the reading produced. And the row needs `&CandidateAlleles`,
which it did not take before.

## 3. The rule, and the part that shipped wrong

**A partial constrains an allele only when every base it showed belongs to a run anchored at a
border.** Four shapes satisfy that:

| the witness | the test |
|---|---|
| one run, flush left | the bases are a **prefix** of the allele |
| one run, flush right | a **suffix** |
| one run, flush both | the read saw every position: **equality** |
| two runs, flush both | the bases divide into a prefix and a suffix, at a point the hole swallowed |

**Everything else is vacuous rather than weak**, and this is where the first version was wrong.
Any other shape leaves a run anchored to neither border, and the bases in that run are
unconstrained — they absorb any disagreement. A three-run witness flush at both borders is
satisfied by taking its prefix and its suffix to be empty; a two-run witness flush at one border,
the same.

**Testing such a witness as though its bases were one contiguous stretch does not weaken the
verdict, it inverts it.** Measured on a six-position locus, a read flush left with a hole at
position 2 that did not reach the right border: the allele agreeing with the read at every
position it actually saw was charged **14 nats**, and the allele disagreeing at a witnessed
position was charged **nothing**. That is the over-restrictive direction spec §5.1 names as the
danger, and it was found by this step's own review, which built the case.

## 4. What a partial is worth

A partial the genotype's carried alleles can all have produced contributes `Σ k_a/P` over them —
a **sum**, since the read could have come from any copy the genotype holds. One none of them can
is charged as an error **with `m = 1`**: a read disagreeing over several positions has no finite
set of wrong outcomes to divide by. The contaminant's half of the mixture is the summed frequency
of the same compatible alleles — a neighbour's DNA shows what this read showed exactly when it
carries one of them.

**A partial compatible with both of a heterozygote's alleles contributes exactly 1**, which is
`ln 1 = 0`: no information, correctly.

## 5. Why an unplaceable read says nothing rather than matching anywhere

A read whose bases cannot be placed could be matched by content — does its sequence occur anywhere
inside the allele? That is sound in one direction only. A carrier's read does occur in its
carrier's sequence, but so may a few bases in an allele it never came from, and short reads in
repetitive sequence make that ordinary rather than rare. **Counting it compatible with everything
costs the information we cannot justify claiming and invents none**: it then contributes the same
to every genotype and cancels, which is §5.3's own account of what a read that saw little should
do.

## 6. What the reviews changed

Two review agents ran on the committed step, each in its own worktree, over reliability (with
mutation testing), errors, naming and smells.

**The inverted verdict of §3 was theirs**, and it is the most serious defect this plan has
produced: a wrong likelihood, not a panic, moving a genotype away from the allele the read came
from by its full charged quality. The rule now decides on the run count as well as the borders.

**And the suite could not tell the implemented rule from the rejected one.** Every fixture spelled
alleles the reference's length and gave every read as many bases as witnessed positions — which is
exactly when the prefix/suffix rule and the positional gather agree. Replacing the two-run arm
with the gather left all 149 tests passing. There is now a fixture where the carrier carries an
**insertion**, so the allele is longer than the reference and the dividing point is not a position
index; it fails on the gather. A second mutation on a different line survived for the same reason —
taking the locus's length from the longest allele rather than from the reference — and has its own
fixture now.

**Three more survivors, all guards nothing exercised.** The `allele.len() < bases.len()` early
return in the split test is load-bearing: without it the two agreements overlap and
`splits_into_a_prefix_and_a_suffix(b"ACTT", b"ACT")` returns `true`, which *adds* copies to a
genotype's score. Both stride assertions on the scratch could be deleted with the module still
green. And the scratch's whole reason to exist — one buffer reused across every sample of a locus
— was never exercised, because all five call sites stood up a fresh one.

**Two gaps in range.** No fixture indexed the copy-share table above two copies, so a partial
compatible with several of a **tetraploid's** alleles — which sums their copies — reached further
into that array than anything tested. And the split helper is a counting form, not a search, so it
is now checked against its own definition over every word of up to four bases against every word
of up to five: 1,953 pairs, exhaustive on a two-letter alphabet.

Five naming findings and the smells findings on `score_partials`' hoist are recorded in
`tmp/review_2026-08-24_calling_likelihood_d1/` and are **not** applied here — they are quality
rather than correctness, and Checkpoint C/D is the place to decide how much of them is worth a
further commit.

## 7. Validation

All in the container:

- `cargo fmt --all -- --check`: clean.
- `cargo clippy --lib --all-features --tests -- -D warnings`: clean. The repo-wide
  `--all-targets --all-features` run is red on `main`, in `examples/ng_duplicated_class_harness.rs`
  and `benches/freebayes_bookkeeping.rs` — pre-existing and out of scope.
- `cargo test`: **4,354 passed, 0 failed, 14 ignored**; **162** in `ng::calling::likelihood`,
  against 141 at C2 and 81 when this milestone run began.

Every mutant named above was re-run against the repaired suite and dies; both files were restored
from byte-identical copies and the restore verified by checksum.
