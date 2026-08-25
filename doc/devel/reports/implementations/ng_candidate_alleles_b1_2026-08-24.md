# ng candidate alleles — B1: the per-sample denominator and the bar

*2026-08-24. Step B1 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md),
branch `ng-candidate-alleles`, on top of `823c7b77`. Design authority:
[`../../ng/spec/candidate_alleles.md`](../../ng/spec/candidate_alleles.md) §1.3, §3, §3.2, §5.1 and
[`../../ng/arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md) §3.1.*

---

## 1. Plan

The first step of Milestone B, and the plan gives it its own commit for a stated reason: **using
the wrong denominator, or forgetting to pool a sample's read-group rows, gives a quietly wrong bar
at every locus and nothing crashes.**

What it builds is the single pass over `CohortObservation::per_sample`. For each covering sample
its *compared reads* are the sum of its support rows — over alleles and over read groups — because
the merge admits only complete observations onto alleles; then each allele it showed, its
read-group rows pooled, is asked whether it reached `max(2 reads, share × compared reads)` against
that denominator. Partial reads are not read.

The pass fills the private `AlleleSummary` that A2 declared: the largest within-sample share, how
many samples cleared the bar, and the cohort read total. Step B2 adds the ranking that reads all
three; nothing in the library calls the pass until C1.

## 2. Assumptions and departures

**The pass fills all three of `AlleleSummary`'s fields, not only the bar's count.** The plan's B2
entry describes the summary as "filled by B1's pass" and then adds `ranks_above`, so the split
between the two steps is *the fold here, the comparison there*. Filling the share and the cohort
total in a second traversal at B2 would be a second pass over the same rows for no gain.

**Read-group rows are pooled by `slice::chunk_by` rather than by
`SampleSupport::pooled_support_for`,** which arch §3.1 names as "the existing method that does
it". The method answers *one* allele, so calling it needs the sample's distinct alleles, and
reaching them means either a scan of the whole table per sample — rows the sample never showed
included — or the same grouping this does. `chunk_by` walks each sample's rows once. It relies on
the merge writing them in ascending `(allele, read group)` order, which `SampleSupport::supported`
documents and `the_rows_are_ordered_by_allele_then_read_group` pins; **the code asserts it rather
than trusting it**, because out-of-order rows would split one allele into two runs and each would
ask the bar with part of the sample's reads.

**A sample whose rows carry no reads at all is skipped.** The merge writes no row for a pair a
sample showed no reads for, so this is unreachable through it; what the arithmetic would otherwise
do is divide by zero and hand B2's ranking a share that is not a number, and arch §2.5 has that
ranking compare with `f64::total_cmp` and no `NaN` branch.

## 3. Changes made

`src/ng/calling/allele_candidates/mod.rs`, two private functions and one import:

- **`compared_reads_of(&SampleSupport) -> u32`** — the denominator, and the only one in the
  module. The sum of the sample's rows over alleles and read groups, the reference's own rows
  included. Its documentation names the three counts the merge keeps that are deliberately *not*
  in it — partials, reads that showed nothing, and reads removed as evidence — and why: none of
  them ever reached an allele, so none can be in a numerator either, and a denominator they
  entered alone would be a bar that rises with a sample's unusable depth.
- **`summarise_alleles(&CohortObservation, MinAltReads, &mut SelectionScratch)`** — the fold.
  Resets the scratch for the locus's table, then walks each covering sample: denominator, then one
  chunk per allele with its read-group rows pooled, updating that allele's cohort read total, its
  best within-sample share, and its count of samples that cleared the bar. The reference is folded
  like any other allele and is asked to pass nothing; whether it is called over is decided
  structurally at C1.

Both carry `#[allow(dead_code)]` with a stated reason — their shipping caller is `select_generic`
at step C1 — for the same reason `AlleleSummary::cleared_the_bar` does: `expect` cannot express
"unused in the library build, used in the test build".

## 4. Tests added

Eight, all on hand-built loci.

| test | what it separates |
|---|---|
| `the_denominator_is_the_samples_compared_reads_and_nothing_else` | the plan's oracle: 9 compared reads counted independently, against a fixture also carrying 6 partial reads, 5 that showed nothing and 3 removed as evidence |
| `one_samples_two_read_groups_are_one_sample` | two lanes of one sample are one sample and one read total, not two of each |
| `the_share_is_one_samples_own_and_the_bar_is_asked_of_each_sample_alone` | the share is a per-sample maximum, and a sample that fails the bar does not enter it |
| `a_sample_showing_only_reference_reads_changes_no_alternatives_summary` | spec §3.2's principle — no term of the bar reads the cohort |
| `a_row_naming_an_allele_the_locus_does_not_hold_is_refused` | a merge bug is refused, not folded into a neighbour |
| `rows_out_of_allele_order_are_refused` | the order the pooling depends on is asserted |
| `folding_a_second_locus_carries_nothing_from_the_first` | the scratch reset, in both directions |
| `a_sample_with_no_reads_produces_no_share_that_is_not_a_number` | the zero-denominator guard |

**The oracle's fixture is built so that every wrong denominator moves the answer.** One sample,
4 reference reads and 5 on one alternative — 3 from one read group and 2 from another — with the
bar at 2 reads or half the sample's compared reads. The bar asks for 5 and the allele has exactly
5. Counting the partials makes it 15 compared reads and a bar of 8; counting the silent reads
makes it 14 and a bar of 7; counting the reads removed as evidence makes it 12 and a bar of 6;
taking the larger read-group row instead of the sum asks 5 of 3. All four fail. Dropping the
reference's own rows makes it 5 compared reads and a bar of 3, which the allele still passes — so
that one error is invisible in the verdict, and the count is asserted directly as well.

## 5. What the review changed

Six agents in six isolated worktrees; full account in
[`../reviews/ng_candidate_alleles_b1_2026-08-24.md`](../reviews/ng_candidate_alleles_b1_2026-08-24.md).

**Two Blockers, and neither was wrong code.** The fold computed the right answer on every input any
agent could build. What was wrong is that four of its properties were pinned by nothing, and eight
mutations passed the suite:

- **the share half of the rule never decided in any fixture**, so deleting the share term outright
  left every test green — and the share is the only half doing any work at 300×;
- **`samples_clearing_the_bar` was asserted as 1 in every place and never above it**, so a fold
  that could not count past one passed — and that count is the cap's deciding key at 3 reads a
  position;
- the share was maximised only in the direction where the largest arrived first;
- nothing pinned what the fold records for the reference, while the doc comment said the opposite
  of what the code did.

**One design conformance failure.** The zero-denominator `continue` is what spec §8 names among the
assertions this step holds in release. Worse than the wording: because it fired before the row
loop, it also stepped over the two conditions the fold *does* assert. It is now an assertion,
separated from the legitimate case — a covering sample whose reads all stopped inside the locus,
which has no rows at all and is stepped over.

Also applied: ascending sample order asserted beside the ascending allele order; the locus region
in all three panic messages; the 157-of-1,707 figure corrected to *more than one* library with its
breakdown; the redundant `dead_code` allow removed and the surviving one documented; the `row` test
helper split so a read group cannot be transposed with a read count; and the parameter renamed from
`bar` to `min_allele_support`, the name the config field carries.

**Four things raised rather than applied**, listed in §4 of the review — the largest being that the
cap's ranking is fed by samples the rule refused, so a one-read sample can decide a truncation.

## 6. Validation

All in the container, on the committed tree:

- `cargo fmt --check` clean;
- `cargo clippy --lib --tests --all-features -- -D warnings` clean;
- `cargo clippy --lib --all-features -- -D warnings` clean, run separately because `dead_code` fires
  there and `--tests` hides it;
- `cargo doc --lib --no-deps` completes, 35 diagnostics all pre-existing and none in this file;
- `cargo test --lib` **4,236 passed, 0 failed, 14 ignored** in 42.5 s, against 4,219 at `823c7b77`.

`cargo clippy --all-targets` is red on `main` with 14 pre-existing errors in five benches and
examples, none in `src/`, which is why the gate is the two `clippy` invocations above.

**Twelve mutations, twelve killed**, each by the test written for it; eight of them survived before
the review's fixes. The table is in §5 of the review report.

## 7. Tradeoffs and follow-ups

- **The fold discards each sample's denominator.** C3 needs it again to ask the rule per
  `(sample, allele)` over the alleles the cap cut. Both versions were built during the review and
  agree to the bit; carrying the denominators would cost one `Vec<u32>` in `SelectionScratch` with
  no reader for two steps, and the cap binds at 23 of 53,935 tomato loci, so the recomputation is
  needed at about one locus in 2,300. **Left to C3 to decide with the measurement in hand.**
- **`examples/ng_candidate_selection_probe.rs` and this fold now disagree by design** — the probe
  asks the rule per `(allele, read group)` row. Its published figures are therefore a lower bound
  on what the shipped rule admits, for the samples carrying more than one read group. D1 deletes
  the probe's copy and is where the two are reconciled.
- **`cleared_the_bar` still carries `#[allow(dead_code)]`.** Its shipping caller remains C1; the
  attribute goes when C1 uses it.
