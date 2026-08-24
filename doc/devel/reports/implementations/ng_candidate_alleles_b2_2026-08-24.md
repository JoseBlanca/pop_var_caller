# ng candidate alleles — B2: the cap's ranking

*2026-08-24. Step B2 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md),
branch `ng-candidate-alleles`, on top of `9aff9c13`. Design authority:
[`../../ng/spec/candidate_alleles.md`](../../ng/spec/candidate_alleles.md) §4, §4.1, §4.2, §8 and
[`../../ng/arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md) §2.5. Milestone B's last
step.*

---

## 1. Plan

A locus is called over at most six sequences including the reference. Above that the list is **cut
to the best six and the locus is still called** — never refused (spec §4.1). B2 builds the
comparison that decides which are the best, over the per-allele fold B1 filled. Four keys, in
order: the largest share of one sample's compared reads the allele took; how many samples cleared
the admission rule for it; the allele's read total over the whole cohort; the bases.

The plan gives it its own commit because **a mis-ordered tie-break is a different truncation at a
minority of loci and nothing fails**. Its oracle: a table built so each level in turn is the one
that decides, plus the same table with its rows reversed, which must not move the answer.

## 2. Assumptions and departures

**Two departures from arch §2.5, both taken on review evidence, both edits to that document that
are not made here.**

- **The name.** Arch calls it `ranks_above`; it is `compare_best_first`. A third-person `-s` verb is
  the shape Rust reserves for `-> bool` (`contains_key`, `starts_with`), and this returns
  `Ordering::Less` for "yes". Two reviewers filed it independently, and a third found the concrete
  cost: `max_by(ranks_above)` returns the **worst**-ranked allele, compiles, and is exactly what
  someone reaching for "the best one" would write.
- **The signature.** Arch takes four positional arguments — two summaries and two base slices. It
  takes two `RankedAlternative` values, each pairing one allele's summary with its own bases. A
  reviewer swapped the two slices at a call site: it compiled, `clippy` was silent, and the test
  still passed, because the bases only decide when all three numeric keys tie. **The mis-pairing is
  invisible at exactly the loci where the ranking does its work**, and the shipping caller is worse
  than the test — step C2 sorts a buffer of table indices, so every argument is an index expression.

**`.then` rather than `.then_with`**, matching the crate's other multi-key comparators; every
argument is an integer or slice comparison, so there is nothing to defer.

**The direction of the fourth key is a free choice.** Spec §4.1 and arch §2.5 both say "then the
bases" without saying which way. Ascending, so the lexicographically smaller sequence ranks above.

## 3. Changes made

`src/ng/calling/allele_candidates/mod.rs`:

- **`RankedAlternative<'bases>`** — one alternative as the ranking reads it: the summary B1 folded
  and the bases that break the last tie, travelling together.
- **`compare_best_first`** — the four keys. The three ranking keys compare right-against-left, which
  is how a descending order is spelled here and matches the crate's other comparators; only the
  bases compare the way they read. Shares use `f64::total_cmp`.
- `SelectionVerdict::Truncated`'s doc comment now names `compare_best_first` instead of forward-
  referencing a name that no longer exists, and `SelectionScratch::ranked_table_indices`' doc gains
  the second order it holds — after the cap has chosen, the survivors are sorted back into the
  merge table's index order, because that is the order they are admitted in.

Both new items carry `#[allow(dead_code)]` naming C2 as their shipping caller.

## 4. Tests added

Twelve, of which one is a property test.

| test | what it separates |
|---|---|
| `the_better_ranked_allele_compares_less` | the direction, which the name cannot carry |
| `the_largest_within_sample_share_decides_first_against_every_other_key` | key 1 against all three others at once |
| `how_many_samples_cleared_the_rule_decides_when_the_shares_tie` | key 2, in the 3-reads-a-position regime where spec §4.1 says it is the only signal |
| `the_cohort_read_total_decides_when_the_share_and_the_sample_count_both_tie` | key 3, production's first key and this ranking's third |
| `the_bases_decide_when_all_three_numbers_tie` | key 4, and that it cannot tie |
| `every_key_decides_a_pair_and_the_row_order_does_not_matter` | the plan's oracle: six alternatives, five adjacent pairs, each key deciding at least one — and the same table reversed |
| `at_three_hundred_reads_a_sample_the_share_alone_decides` | the deep end of the committed range |
| `at_one_sample_the_share_ranking_and_the_cohort_total_ranking_agree` | the thin end of the cohort axis |
| `a_shallow_homozygote_outranks_a_deep_heterozygote_on_the_first_key` | a cohort of mixed depth, which spec §4.1's argument does not cover |
| `a_table_above_the_cap_ranks_its_survivors_into_the_front` | seven alternatives, the first width at which the cap bites |
| `two_shares_equal_to_the_bit_fall_through_to_the_next_key` | equal shares, and an honest note that no fixture separates `total_cmp` from `partial_cmp` |
| `the_ranking_is_a_strict_weak_ordering` (property) | asymmetry and transitivity, which `sort_unstable_by` requires and does not check |

**Every pairwise fixture sets the bases against the expected answer**, which is a fix rather than a
style: see §5.

## 5. What the review changed

Four agents in four isolated worktrees; full account in
[`../reviews/ng_candidate_alleles_b2_2026-08-24.md`](../reviews/ng_candidate_alleles_b2_2026-08-24.md).

**The step's own oracle did not test what its doc comment claimed.** With the second and third keys
swapped in the comparator, `every_key_decides_a_pair_and_the_row_order_does_not_matter` still
passed: the pair whose deciding level was meant to be the sample count had *equal* cohort totals, so
the third key was a no-op there and a comparator consulting it first reached the same answer. One
number in the fixture fixes it, and the mutation now fails the oracle.

**And the base tie-break silently agreed with every numeric key.** In all seven pairwise fixtures
the expected winner's bases also sorted first, so replacing the entire function body with
`left.bases.cmp(right.bases)` — no numeric keys at all — passed seven of the eight tests. Every
fixture now sets the bases against its answer; that mutation fails nine.

**Four claims in the new prose were wrong, and every one was about my own fixture.** The 300-read
test said both lower keys were set against the answer where the cohort total was set with it; the
six-row table's reading line said an allele "loses on the bases" where it wins; the first-key test
said the share "is the only key that separates them" where all three do; and the `total_cmp` test
attributed a `None` return to `partial_cmp` on two equal shares, which returns `Some(Equal)`. All
four corrected, and the fixtures changed where the prose was describing what it should have said.

Also: the rename and the two-argument signature above; a redundant `assert_ne!` implied by the
`assert_eq!` three lines above it; and three tests added for coverage the review named — one sample,
mixed depth, and a table wider than the cap.

**Three things raised rather than applied**, in §4 of the review report. The largest is what the
already-open key-1 question costs *through* the ranking: because key 1 is compared before keys 2 and
3 can speak, one read from a sample that failed the admission rule can no-call 40 samples that
passed it.

## 6. Validation

All in the container, on the committed tree:

- `cargo fmt --check` clean;
- `cargo clippy --lib --tests --all-features -- -D warnings` clean;
- `cargo clippy --lib --all-features -- -D warnings` clean;
- `cargo test --lib` **4,248 passed, 0 failed, 14 ignored** in 45.8 s, against 4,236 at `9aff9c13`.

**Eight mutations, eight killed.** Before the review's fixes, two of them were survivors — the
bases-only comparator, and the swapped key order against the oracle. The table is in §5 of the
review report.

## 7. Tradeoffs and follow-ups

- **Three `#[allow(dead_code)]` attributes now stand in the module** — `compare_best_first`,
  `summarise_alleles` and `AlleleSummary::cleared_the_bar` — and C2 removes all three at once. A
  reviewer's C1–C3 prototype needed exactly that.
- **`SelectionVerdict::Truncated { dropped: u16 }` cannot count what the cap can cut.** A reviewer
  demonstrated the conversion failing at 70,001 alternatives on one locus. It needs 65,536 distinct
  sequences at one position, which neither benchmark approaches, but spec §4.1's whole argument is
  "truncate, never refuse" and this is the one input on which C2 would refuse by panicking. Raised
  at Checkpoint B; it is an edit to arch §2.2.
- **The surviving alternatives reach the VCF in the merge table's order, not in rank order.** The
  ranking decides *which* survive; admission is in table order (arch §3.1). True and nowhere stated,
  now recorded on `ranked_table_indices`.
