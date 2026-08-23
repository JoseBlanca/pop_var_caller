# ng genotype prior — the fit now reports how far its answer is, instead of claiming it matched

*Fix-forward report, 2026-08-23. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Repairs `0b019e0d`, which shipped before the review fan-out
that every plan step gets; the owner folded that commit into Milestone E's review, and this is
what it found.*

## 1. The defect

`SeedRegime::FittedSpectrum` carries a marker saying whether the two starting numbers the fit
reads off the panel's allele-frequency spectrum reproduce what was measured. **Nothing compared
them.** The marker was set to *reproduced* whenever two unrelated conditions held: the search
finished in the interior of its range, and no allele-count class came back predicted at exactly
zero. Neither says anything about how close the answer is.

Measured, on panels the fit will meet:

| measured spectrum | fitted pair's prediction | old marker |
|---|---|---|
| 5 classes, `[0.05, 0.45, 0, 0.45, 0.05]`, 2 individuals | `[0.099, 0.244, 0.313, 0.244, 0.099]` | *reproduced* |
| 53 classes, all weight at two middling frequencies, 26 individuals | shares **4 parts in 100** of its mass with the measurement | *reproduced* |

Both are the shape spec §4.1 names as the one a two-parameter family cannot hold, so both are
exactly what the marker existed to catch. The commit message, the type's own documentation and
`arch/calling_priors.md` all claimed it caught them.

## 2. The fix, and why it reports a distance rather than a verdict

The marker is now a distance: **the Kullback–Leibler divergence of the measured spectrum from the
fitted pair's prediction, in nats**, with zero meaning the family reproduced the measurement
exactly. It names no threshold, and that is the point — **nobody has measured how far off the
starting pair has to be before a genotype moves**, so classifying was the part that had to go, not
the checking. Whoever reads a run's output decides what is too far, and the number is there to
decide with.

**It costs nothing.** The fit's objective is already the measured spectrum's own entropy minus
this divergence — the type's own documentation said so before the fix, and a test
(`the_winning_score_is_the_spectrums_own_entropy`) already exercised the identity. So subtracting
the winning score from that entropy gives the distance with no prediction at all. The old check
*did* predict once more, so **a fit is back to the 399 predictions it cost before that commit**,
where the commit message quoted 400.

Whether the search ran out of range is carried separately, because it is not derivable from the
distance: a pair pinned against a bound can still predict the measurement well, and the invariant
cohort in the tests does exactly that — at the search limit, at a divergence of 1.0e-9 nats.

**Reference values, all measured in this module's tests:**

| spectrum | divergence |
|---|---|
| any shape the family can hold, over three panel sizes and three inbreeding coefficients | worst **1.1e-9 nats** — what the search's 1% resolution leaves behind |
| the five-class bimodal panel above | **0.481 nats** |
| the 26-individual bimodal panel above | **3.153 nats** |
| a panel at inbreeding exactly 1 holding heterozygotes | **above 10 nats** — the objective charges the impossible classes `ln(PROBABILITY_FLOOR)` |

The last row is what the old `Unreproducible` variant caught, and it is worth noting what the
distance does better: that variant fired on **any** class with weight above zero that the pair
could not produce, so a class carrying 1 part in 10,000 flipped the whole run's marker. A
divergence weights the class by how much of the panel actually sits in it, which is the right
answer to the same question.

## 3. Four further findings against the same commit, applied here

- **`fill_locus_concentration` was not `#[must_use]`**, so a loop could discard the checked
  concentration it returns and go on to use its own raw buffer — the hole the commit was written
  to close, left half open. Adding the attribute immediately caught **six call
  sites in the commit's own tests** that discarded the value, three of them in tests that assert a
  panic and so never reach the return at all.
- **An inline comment said the SNP-versus-indel split belongs to the projection**, which is the
  opposite of the decision the same commit records forty lines above and in the spec.
- **`ProjectionFit::predictions` was documented as 399** while the commit had made it 400. It is
  399 again, and the doc now says why.
- **`an_empty_seed_is_refused` no longer tested what it was named for.** Since the seed arrives as
  a checked type, `Concentration::new(&[])` panics while the argument is evaluated, so
  `fill_sample_concentration` is never entered; the test had become a duplicate of
  `an_empty_concentration_is_refused` under a name describing a call it did not make. Deleted, and
  its reasoning moved onto the test that does cover the behaviour.

## 4. What the tests pin

- `a_spectrum_the_family_can_hold_scores_at_effectively_zero_divergence` — the baseline, worst
  1.1e-9 nats over three panels, and that none of them is at the search limit.
- `a_spectrum_the_family_cannot_hold_scores_far_from_it` — **the test the old marker could not
  pass**: both bimodal panels above, at 0.481 and 3.153 nats.
- `a_spectrum_no_pair_can_produce_is_marked_rather_than_answered` — above 10 nats at inbreeding 1.
- `a_fit_that_stopped_at_the_edge_of_its_range_says_so` — at the limit **and** at 1.0e-9 nats,
  which is why the two facts are carried separately.

## 5. Gates

Green in the container: `cargo fmt --check`, `cargo clippy --lib --tests --all-features -D
warnings`, `cargo test --lib`.
