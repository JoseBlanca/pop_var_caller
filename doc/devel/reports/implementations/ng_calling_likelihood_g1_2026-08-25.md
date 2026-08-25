# ng read likelihoods — G1: the censored term, and a property the specification does not have

*Implementation report, 2026-08-25. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step G1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md). **This completes
Milestone G.***

## 1. What it is

A read that entered a repeat tract and **ran off its own end** does not say what length the sample
carries there. It says the tract is **at least** as long as the stretch it got through — censored
data, in the statistician's sense. `StutterSubstitutionEmission::censored_emission` puts a number
on that, and until this step it was `unimplemented!()`.

The number is the letter match over what the read witnessed, times the probability that the tract
came out at least that long (spec §5.2). The second factor is new and is where all the work is.

## 2. Two forms, and which candidate gets which

**On a pure tract — whole copies of the motif end to end — the letters factor out exactly.**
Stretching a pure tract appends or trims at its end, so every stretching agrees on the first `ℓ`
bases and the substitution term comes out of the sum. What is left is a tail:

```text
Lr(a read that got through ℓ bases | candidate a)
    =  subst( the ℓ bases | first ℓ bases of a's own tiling )  ×  P( length ≥ ℓ | a )
```

**On an interrupted tract they do not.** The stretchings put the interruption at different
offsets, so the first `ℓ` bases genuinely differ between them and the code sums change by change,
slicing each realisation to the witnessed length. Spec §5.2 bounds what factorising here would
cost at `log(3(1−ε)/ε)` per distinguishing base — 6.4 nats at an error rate of 1 in 200 — and says
to pay for the sum instead. `is_a_pure_tract` is the one-line predicate that chooses.

## 3. Where the tail lives, and why it is not in `calling/`

`P(length ≥ ℓ | a)` is summed by two new methods on `StutterModel`
([`alignment/stutter.rs`](../../../../src/ng/alignment/stutter.rs)):

- **`reachable_length_changes`** — an iterator over every length change a read of this candidate
  could actually show, each paired with its probability. At most 41 of them, whatever the tract's
  length, because the two cutoffs bound the support.
- **`probability_at_least_this_much_longer`** — that support, filtered by one inequality and
  summed.

**They are there and not under `calling::likelihood` because the map from a part-repeat step back
to a length change is the inverse of the one `StutterModel::probability` applies**, twenty lines
above it. Two inverse functions in two modules have nothing holding them together; here one test
round-trips them.

**The reachability rule is the second reason and the weaker one, because it is genuinely stated
more than once.** A read of a candidate must still show a repeat, so at most `repeat_count − 1`
repeats may be contracted away, and at most `(repeat_count − 1) · (period − 1)` re-indexed steps
on the part-repeat branch. `unreachable_mass` states that boundary from the other side, and the
two now share one function, `contractable_repeats`, so they cannot drift. A **third** statement
lives in `enumerate_placements` at a different grain — it sees the tract's runs where the
distribution sees only the total — and that gap is the recorded open question rather than a copy
to delete. The review's own words for it: the move would put the distribution-grain statement in
two files and leave the run-grain one where it is, so the count goes from two to three rather than
staying at two.

**The three are now tied by tests, which they were not.** Before this step, moving the
distribution-side bound broke no test in `calling::likelihood`, and moving the byte-level rule
broke exactly one test that pins it with hand-written literals. Now
`the_placements_and_the_support_agree_on_a_pure_tract` fails on the byte-level mutation and
`the_handed_out_support_totals_what_the_loss_reports` on the distribution-side one; moving
`contractable_repeats` alone fails 11 tests across both files.

**An agent ran the counterfactual.** It moved both methods into `calling::likelihood` and
confirmed the result compiles and passes — at the cost of 24 duplicated lines of test oracle,
because the base-pair walk is also `unreachable_mass`'s oracle and would have to stay behind.

## 4. Two departures from the plan, and the reason for each

**The tail is a term-by-term sum, not the telescoped closed form the plan names.** Spec §12's
twelfth test asks that where the constraint admits exactly one length change, the censored
likelihood equal the complete likelihood at that change **bit for bit** — and calls that "the test
of the tail arithmetic rather than of the tolerance". A tail summed from `StutterModel::probability`'s
own terms satisfies that by construction; a telescoped `(1 − g)^(a−1) − (1 − g)^b` is equal only
algebraically and fails it. The sum is still exact and finite — the cutoffs bound the support, so
nothing is truncated by taking it — and it is at most 41 terms. The plan's word was *closed-form*;
the property the specification actually requires decides against it.

**Spec §12's thirteenth test is not implemented as stated, because as stated it is false.** §5.2
says *"A partial is always less discriminating than a complete observation of the same bases,
because a tail probability varies less between candidates than a point probability does"*, and §12
asks for that without restriction. Section 5 below is the counterexample, with its sizes. Two tests
carry the honest form.

## 5. The finding: a partial can out-discriminate a complete read, and by a lot

**Where the two candidates straddle the stretch the read got through, the partial separates them
further than the complete read does.** Measured, at motif `CA`, an error rate of 1 in 1,000, and
the contraction-biased fitted row the module's tests use (slippage level 2 in 100, 83 in 100 of it
contraction, fall-off 0.35):

| | separation between a 4-repeat and a 6-repeat candidate |
|---|---|
| a **complete** read of 10 bases | **1.586 nats** |
| a **partial** read that got through 10 bases | **5.661 nats** |

**Why, in one sentence each.** The complete read needs a one-repeat expansion under the 4-repeat
candidate and a one-repeat contraction under the 6-repeat one, so what separates them is the
direction ratio alone. The partial read needs an expansion under the 4-repeat candidate and
**nothing at all** under the 6-repeat one — a 12-base tract is already at least 10 — so it collects
the same-length share, which is most of the distribution.

**That is real information, not a defect.** A lower bound rules out everything below it, which is
the evidence spec §5.1 turned these reads on to collect in the first place. What is wrong is the
claim about it.

**Checked independently, and the check named what §5.2 was reaching for.** A review agent
reproduced both numbers from an oracle that walks base pairs calling `StutterModel::probability`
directly — touching neither new method — and agreed with the code to better than one part in
10¹². Its diagnosis: the theorem behind §5.2's sentence is the data-processing inequality, which
bounds **expected** discrimination. §5.2 states the pointwise version of an averaged result, and
the pointwise version is false. **Nothing here depends on the parameters.**

**The safety property §5.1 actually needs is a different one, and it holds.** A censored read
never scores *below* the complete read of the same bases, so it can never be mistaken for
evidence of a short allele — which is the trap §5.1 names. That is what
`a_censored_read_is_at_least_as_likely_as_the_complete_read_of_the_same_bases` pins.

**The restricted form is very nearly true and is what the test now pins.** Among candidates the
read **outgrew** — both shorter than what it saw — the geometric's memorylessness makes the tail
proportional to the point probability and the two separations come out almost equal. Over a grid of
two parameter rows, three read lengths and six candidate pairs they differ by at most **0.043
nats** against a separation of **3.149 nats** — 1.4 parts in a hundred — and the partial is the
larger at that cell, so even here "no larger for the partial" is false. What breaks the exact
equality is the part-repeat branch, whose re-indexing means its terms do not line up with the
whole-repeat ones.

**Both facts are asserted by tests rather than stated here**
(`a_censored_read_out_discriminates_a_complete_one_where_the_candidates_straddle_it` and
`a_censored_and_a_complete_read_separate_outgrown_candidates_alike`), so they cannot go stale
silently.

**What the owner has to decide, and it does not block H1 or H2.** Spec §5.2's sentence and §12's
thirteenth test are wrong and should be corrected — the recommendation is to replace the "always
less discriminating" claim with the two properties that are true: a censored read never scores
below the complete read of the same bases, and among candidates the read outgrew the two separate
them by the same amount to within 1.4 parts in a hundred. **Nothing in the model changes either
way**: the formula §5.2 specifies is implemented exactly as written, and only the claim about its
behaviour moves. Left as a stop-and-ask because editing the design documents is not this loop's
to do.

## 6. What the tests pin

**In `alignment/stutter.rs`, five:**

- `the_handed_out_support_totals_what_the_loss_reports` — the iterator's total against the
  independent base-pair walk and against `unreachable_mass`, over 4 parameter rows × 6 periods × 6
  tract lengths, plus a no-duplicates check.
- `the_censored_tail_and_its_complement_make_the_reachable_total` — **spec §12 test 12**, at
  **8,364 floors**. The tail comes from the model's own support; the complement from the base-pair
  walk, so the two halves are not two readings of one expression. The shortfall from 1 that the
  sweep spans is asserted rather than described: **widest 2.0618 parts in 100** — a mononucleotide
  tract of one repeat at a one-step share of 0.01, which can neither contract nor reach the
  part-repeat branch — against the 1e-12 the identity is checked to, so a version of this test that
  compared against 1 fails rather than passing by luck. **Narrowest exactly zero**, which is
  `unreachable_mass`'s clamp rather than a perfect total.
- `a_tail_at_either_end_is_the_whole_distribution_or_nothing` — the two ends, which the sweep
  averages over.
- `a_tail_of_one_change_is_the_probability_of_that_change` — the bit-for-bit reduction, on the
  distribution alone.
- `part_repeat_bp_diff_inverts_the_re_indexing_the_distribution_applies` — the step-to-base-pair
  map against the one `probability` applies, so a broken re-indexing reports itself at the line
  that is wrong rather than as a distribution total that disagrees.

**In `calling/likelihood/ssr_emission.rs`, twelve** — five written with the step and seven added
by the review:

- `a_censored_read_of_one_admissible_length_scores_what_the_complete_read_does` — **spec §12 test
  12's second half**, `to_bits()` equality, on three pure tracts and one interrupted one, so both
  branches are covered.
- `a_censored_read_is_at_least_as_likely_as_the_complete_read_of_the_same_bases` — the cheapest
  check that the tail sums at-or-above rather than something narrower.
- `a_censored_and_a_complete_read_separate_outgrown_candidates_alike` — **spec §12 test 13**, in
  the form that holds.
- `a_censored_read_out_discriminates_a_complete_one_where_the_candidates_straddle_it` — the
  counterexample, with its two sizes.
- `a_censored_read_past_every_stretching_scores_nothing` — falls to the row's outlier term, as an
  impossible complete length does.

## 7. What the review found

**Seven category agents, each in its own worktree.** They ran 33 mutations between them and
returned **two Blockers, seven Majors** and a long tail of Minors. Everything below is fixed in
this commit unless it says otherwise.

### The two Blockers were both tests that could not fail

- **`is_a_pure_tract`'s length check was untested.** Deleting the `is_multiple_of` conjunct
  left the whole suite green, and the mutant is not a no-op: `chunks_exact` drops the
  remainder, so `CAGCAGTT` under motif `CAG` comes back *pure* when it plainly is not. Measured
  on that tract, a censored read of `CAGCAGCAG` scores **3.369e-3** down the interrupted route
  against **3.939e-10** down the pure one — a factor of eight million, arriving with no panic.
  Nothing in the file reached it: the one non-multiple candidate that gets this far fails
  *both* conjuncts, so the length check was never the one deciding.
- **Both prefix cuts on the interrupted branch could be deleted with the suite green.** The
  exact-sum branch was reached by exactly one fixture, and in it the read was exactly as long
  as the single stretching admitted — so both cuts were full-length no-ops. The branch's three
  distinctive behaviours — cutting a proper prefix, the part-repeat arm, and summing more than
  one term — had no test at all.

### The Majors

- **A doc quoted a loss range its own cited test refutes.** The tail's documentation said the
  complement identity leaves "2 parts in a million to 2 in a thousand". Those figures are real
  but belong to one cell — `unreachable_mass` introduces them for the contraction-truncation
  term alone at four hexamer repeats. The test cited beside them measures **0 to 2.06 parts in
  a hundred**, ten times the stated ceiling at one end and exactly zero at the other, where the
  same sentence's "wrong by orders more than any tolerance" is simply false.
- **The trait's own documentation still carried the claim section 5 disproves**, five lines
  above the method that disproves it.
- **Nothing pinned that a pure candidate takes the factorised route.** Forcing every candidate
  down the exact sum passed the whole suite, though the two routes differ in the last bits on
  241 of the cells swept — so a bitwise assertion separates them where a tolerance cannot.
- **The interrupted branch re-split the tract once per length change.** Measured at 21
  segmentations per (observation, candidate) pair where `emission` does one, and 15.6–21.1 µs
  against the pure branch's 141–193 ns, of which 92% was the re-splitting. The scratch type's
  own doc sentence — "nothing allocates per observation per candidate" — was false as written.
- **A test named a failure it cannot catch.** `a_censored_read_is_at_least_as_likely_...` said
  it catches a tail that dropped its own floor. It does not: dropping the floor makes the tail
  *larger*, which an inequality in that direction accepts by construction. Replacing the floor
  with `i64::MIN` leaves that test green and fails four others.
- **Every censored fixture was an exact tiling prefix**, so the letter half of the product was
  at its maximum in every censored call in the suite and never varied.
- **`emission` and the support disagree about what a tract can reach.** See section 9 — this one
  is *not* fixed here.

### What changed structurally, and why it was worth it

The placement dispatch was written twice, and the two copies had already drifted — the
unreachable-case comment existed in two wordings saying different things. It is now one
function, `letters_over`, and the review's own question ("does extracting obscure the
difference between cutting to the witnessed length and not cutting?") has the answer **no,
because that difference is not one**: in `emission` the length change is derived from the
observation, so cutting to the observation's length is the identity there. The complete read is
the degenerate case of the censored one. `censored_emission` went from 90 lines to 55, and the
candidate's runs are now found once per call rather than once per length change.

**Five mutations re-run against the repaired tests, all now killed:** dropping the length check
(1 failure), forcing the exact sum (1), removing the prefix cut (2), letting the byte-level rule
empty the tract (2 — including the new test that ties it to the distribution's rule, which
nothing did before), and moving the distribution-side contraction bound (11, now spanning both
files where it previously touched no likelihood test at all).

### Two judgements the review returned that I did not act on

- **The module placement was challenged and upheld.** An agent moved both methods into
  `calling::likelihood` and confirmed it compiles and passes — at the cost of 24 duplicated
  lines of test oracle and, the part nobody had named, separating the step-to-base-pair map
  from its own inverse inside `probability`. The documentation now leads with that argument
  instead of the weaker one it had.
- **The tail's floor prunes no work.** An agent wrote the obvious optimisation, confirmed it
  bit-identical, and measured it **slower** on the common case. Recorded so it is not tried
  again.

## 8. Validation

Run in the dev container on this worktree:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features --tests -- -D warnings` — clean. (`--all-targets` stays red on
  `main` for unrelated reasons in `examples/ng_duplicated_class_harness.rs` and
  `benches/freebayes_bookkeeping.rs`.)
- `cargo test --lib` — **4,410 passed, 0 failed, 14 ignored**, against a baseline of 4,392;
  **195 in `ng::calling::likelihood`** against 182, and **32 in `ng::alignment::stutter`**
  against 27.

## 9. What is still open, unchanged by this step

- **`emission` scores 22 length changes the support calls unreachable, and repairing it is a
  decision about the model.** All 22 run one way — `emission` places mass the support refuses,
  never the reverse — and every one is a part-repeat contraction that would leave the tract
  below one repeat. A two-base `CA` tract scores a one-base read at 5.390e-4 while
  `unreachable_mass` counts that same mass as unplaced, so the number the row uses to keep
  candidates comparable is describing a different function from the one doing the scoring.
  The cause is that `emission`'s part-repeat branch asks only whether the distribution gives
  the change a non-zero probability, and that question does not know the tract's length.
  **It is milestone F's code, not this step's**, and the repair moves scores at one- and
  two-repeat tracts, so it is left as it is and pinned by
  `the_scoring_and_the_support_disagree_only_where_the_open_question_says_they_do`, which
  asserts the count and the direction.
- **`unreachable_mass` understates the loss on interrupted tracts** (Checkpoint F handoff, item 3).
  The censored term inherits it exactly: its exact-sum branch takes `enumerate_placements`' answer,
  which sees the runs, while the tail takes the distribution's, which sees only the total. The two
  agree on every pure tract and can differ on an interrupted one. Papering over it here would have
  made the two branches disagree with `emission`, which is worse.
- **The oracle's `censored_emission` is still `unimplemented!()`**, with its existing message: Model
  B scores complete observations only. Its independence from Model A is in how it explains a read's
  *length*, and a censored read's length is not explained but bounded, so a censored Model B would
  have to truncate to the witnessed prefix exactly as Model A does and would stop being a second
  opinion. Plan step G1 asks for the seam's method; F3's cross-model check is scoped to complete
  observations.
