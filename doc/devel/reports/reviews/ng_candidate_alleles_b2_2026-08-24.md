# ng candidate alleles — B2: review and the fixes applied

*2026-08-24. Step B2 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md).
Reviewed at `9aff9c13` plus the step's working-tree diff, four agents in four isolated worktrees.
Implementation report: [`../implementations/ng_candidate_alleles_b2_2026-08-24.md`](../implementations/ng_candidate_alleles_b2_2026-08-24.md).*

---

## 1. Which categories ran

`reliability`, `naming`, `idiomatic`, and **design fidelity** — the fourth again asked whether the
steps that follow could be built on this one, and told to answer by writing them.

`errors` was not dispatched: the step adds a pure comparison with no fallible path, no `Result`, no
panic and no I/O. `smells` was not dispatched: one function and one two-field struct, and the
duplication question against the measurement harness was already settled at B1.

Each agent was told what the author had already mutated, so none spent a round re-finding it.

## 2. What was actually wrong

### 2.1 The step's own oracle did not test the key order it claimed to — Major

The plan gives B2 its own commit because "a mis-ordered tie-break is a different truncation at a
minority of loci and nothing fails", and
`every_tie_break_level_decides_one_pair_and_the_row_order_does_not_matter` was written to catch
exactly that. **With the second and third keys swapped in the comparator, it still passed** — only
the single-purpose sample-count test failed.

The cause was one number. The pair whose deciding level was meant to be the sample count was
`(0.5, 3, 8)` against `(0.5, 2, 8)`: their cohort totals were **equal**, so the third key was a
no-op on the only pair the second key was there to decide, and a comparator consulting the totals
first reached the same answer anyway.

**Fixed** by setting that allele's cohort total against the expected answer, 8 → 30. The expected
list is unchanged, and the swapped-key mutation now fails the oracle as well as its sibling. The
doc comment gained the sentence that makes the number legible: "`AG` carries the larger cohort
total, so a ranking that consulted the totals before the sample count would swap them."

### 2.2 The base tie-break silently agreed with every numeric key — Major

Across all seven pairwise fixtures the expected winner's bases also sorted first — always `C` or
`AC` beating `G` or `AG`. So **replacing the whole function body with `left.bases.cmp(right.bases)`,
keeping no numeric key at all, passed seven of the eight tests.** All the discriminating power for
three separate rules sat in the single six-row table, which is a single point of failure for them.

**Fixed** by setting the bases against the expected answer in every pairwise fixture, and saying so
where a reader will find it. That mutation now fails nine tests. The exception is the test about
the bases themselves, where they have to decide.

### 2.3 Four claims in the new prose were wrong — Major in aggregate

Every one was a statement about my own fixture, which is the defect class the implementation skill
names as the most reliable this loop produces. The naming agent checked each at its source:

- the 300-read test said all the lower keys were set *against* the answer; the cohort total, 151
  against 48, was set **with** it, so a share-free ranking would have passed;
- the six-row table's reading line said `CA` "loses on the bases" — it **wins**, and the expected
  list ranks it above `CC`;
- the first-key test said the share "is the only key that separates them" — the sample count and
  the cohort total separate them too, which is that paragraph's own next point;
- the `total_cmp` test attributed a `None` return to `partial_cmp` on two equal shares, which
  returns `Some(Equal)`; the fixture would pass unchanged against `partial_cmp`.

**All four corrected**, and the fixtures changed where the prose was describing what it should have
said: the 300-read fixture's error-level allele now carries 250 cohort reads against the
heterozygote's 151, so all three lower keys really do oppose.

The same agent verified the figures that *were* right rather than only reporting the wrong ones:
the 15-against-240 and 0.5-against-0.01 scale argument closes arithmetically and matches spec §4.1;
"two thirds at 3 reads" matches spec's "near 0.67"; and all three external code citations
(`enforce_max_alleles` really is a stable sort on `Reverse(cohort_count)`, `bam2bcf.c:975`,
`HaplotypeGenerator.cpp:183`) check out at source.

### 2.4 Four positional arguments could be mis-paired invisibly — Major

`ranks_above(&AlleleSummary, &[u8], &AlleleSummary, &[u8])` pairs each allele's summary with its own
bases **only by position**. The idiomatic agent swapped the two slices at a call site: `clippy
--lib --tests -- -D warnings` was clean and the test still passed, because the bases only decide
when all three numeric keys tie. **The mis-pairing is invisible at exactly the loci where the
ranking does its work.** The shipping caller is worse than the test — the design-fidelity agent's
C2 sorts a buffer of table indices, so every argument is an index expression with `left` and `right`
interleaved.

**Fixed** by `RankedAlternative { summary, bases }`, two arguments instead of four. The agent wrote
it, formatted it and ran the gate before proposing it: `clippy --lib` and `--lib --tests` clean,
49 tests passing. It needs no `allow` of its own, because the `allow` on the comparator makes it a
live root.

*This changes the signature arch §2.5 prescribes, and that document is not edited here.*

### 2.5 The name asks a question and answers `Less` for yes — Minor, filed by three agents

`ranks_above` has the shape Rust reserves for `-> bool` predicates, and returns `Ordering`. The
patch supplied its own evidence: the doc comment spent a whole paragraph undoing the name, and a
test existed only to pin the direction, its first line reading "which the name alone does not give".

The design-fidelity agent found the concrete cost while writing C2. The sort is fine —
`sort_unstable_by(ranks_above)` gives best-first and reads correctly. The trap is the neighbouring
idiom: **`max_by(ranks_above)` returns the worst-ranked allele**, compiles, and is what somebody
reaching for "the best one" would write. Its failure mode is keeping the wrong five alleles at a
truncated locus — silent, and rare enough (23 of 53,935 tomato loci) not to show up in a smoke test.

**Renamed to `compare_best_first`**, which states the direction in the name. Cost, all paid here:
one committed doc comment in the module, six test call sites. Arch §2.5 and the plan still say
`ranks_above` in prose; that is an owner edit, flagged below.

### 2.6 Smaller things, all applied

- **A test asserted a `partial_cmp` guard it did not provide and no fixture could.** Substituting
  `partial_cmp(..).unwrap()` leaves every test in the module green, because a `NaN` share cannot
  reach the comparator — B1's fold asserts a non-zero denominator for every sample with rows. The
  test now says that outright: `total_cmp` is the right spelling because it cannot panic if that
  assertion is ever relaxed, but its value here is a guarantee, not a behaviour any fixture sees.
- **`assert_ne!(.., Ordering::Equal)` was implied** by the `assert_eq!(.., Ordering::Less)` three
  lines above it on the identical call. Removed.
- **`use std::cmp::Ordering` sat below the `crate::` block** where sibling `ng` modules put `std`
  first. Moved.
- **`.then_with` over integer comparisons** where the crate's other multi-key comparators use
  `.then`. Changed; there is nothing to defer.
- **Three coverage gaps named by two agents, all now tested**: one sample (where this ranking and
  production's provably coincide, which bounds what spec §4.1's argument buys); a cohort of mixed
  depth; and a table of seven alternatives, the first width at which the cap bites — every other
  fixture was six rows or fewer.
- **A property test for the strict weak ordering**, which `sort_unstable_by` requires and does not
  check: given a comparator that is not one, the standard library may leave the slice in any order,
  and spec §8 requires byte-identical output at any worker count. `proptest` was already a
  dev-dependency.
- **`ranked_table_indices`' doc described only its first state.** Writing C2 showed it necessarily
  has two: ordered by the ranking while the cap chooses, then sorted back into merge-table order,
  because admission is in table order. Extended.

## 3. Checked and found sound

- **Fidelity: the four keys are spec §4.1's, in that order, with that meaning.** The design-fidelity
  agent found no place where code and documents disagree about them, and verified the two supporting
  claims in the doc comment rather than assuming them: no two entries can compare `Equal`, because
  the merge's `AlleleTable` keys its allele vector by the bytes themselves; and no share can be
  `NaN`, because B1's fold asserts a non-zero denominator.
- **The comparison is a valid strict weak ordering.** The reliability agent checked antisymmetry and
  transitivity exhaustively over a set including `NaN` and `-0.0` shares and a prefix pair.
- **Nothing about B2's shape forces a change to C1, C2 or C3.** The design-fidelity agent wrote all
  three in `generic.rs`, the file arch §3.1 puts them in, plus 13 tests including the plan's C1
  round-trip oracle and its second C3 oracle. It built first try and ran 62 tests green.
- **The `ranked_table_indices` buffer works without a per-locus allocation**, which was the one
  question that could have forced a shape change. Destructuring `&mut SelectionScratch` into its two
  private fields gives disjoint borrows, so C2 sorts the index buffer while reading both the
  summaries and the bases. `select_generic` must live inside the `allele_candidates` module tree for
  that, which arch §3.1 already requires.
- **`f64::total_cmp`, and right-against-left for a descending order**, both match the crate's
  existing comparators.

## 4. Raised, not applied — decisions for the owner

Put to the owner at Checkpoint B. None blocks the commit.

1. **What the open key-1 question costs through the ranking, which is more than it looked.** B1's
   fold maximises the within-sample share over *every* sample, including ones that did not clear the
   admission rule. Because key 1 is compared before keys 2 and 3 can say anything, **a share of 1.0
   contributed by a single-read sample is lexicographically dominant over the whole rest of the
   ranking** — key 2, the only signal at 3 reads a position, cannot overturn it. The design-fidelity
   agent built the smallest case and it passes: an allele 40 samples carry at 150 of 300 reads each,
   6,000 cohort reads, is cut in favour of an allele one sample carries — because a second sample's
   single read landed on it. C3's second count then emits **all 40 carriers as missing**. Two
   candidate fixes, both one line inside B1's existing loop: maximise over samples that cleared the
   rule for that allele, or over samples whose compared reads reach the rule's floor. Either way
   spec §4.1's sentence should say which samples the maximum is over.
2. **A cohort of mixed depth, which spec §4.1's range argument does not cover.** Its two halves each
   describe a cohort at one depth. Where the regimes meet, a homozygote scores 1.0 whatever its
   depth and a heterozygote about 0.5 whatever its depth, so a 3-read sample's allele — three reads
   in the whole cohort — outranks every 300-read heterozygote and is the last thing a binding cap
   cuts. Not obviously wrong (three agreeing reads *are* evidence), but a reader of §4.1 would
   conclude the opposite. Recorded as a test rather than argued about; the paragraph is spec §4.1's
   and §7's to gain.
3. **`SelectionVerdict::Truncated { dropped: u16 }` cannot count what the cap can cut.** Demonstrated
   with a locus of 70,001 alternatives, where the conversion panics. Nothing in the merge bounds a
   locus's allele table at 65,536 — `CandidateAlleles::admit`'s refusal at that width guards the
   *candidate* table, which the cap holds at six. It needs 65,536 distinct sequences at one position
   and neither benchmark comes close, so it is small; it is raised because spec §4.1's whole argument
   is "truncate, never refuse" and this is the one input on which C2 would refuse, by panicking.
   Widening to `u32` is the honest fix and costs two bytes per locus verdict; it edits arch §2.2.
4. **Arch §2.5 and the plan still name `ranks_above` and the four-argument signature.** Both changed
   here on review evidence (§2.4, §2.5 above), and this loop does not edit design documents.

## 5. Validation after the fixes

All in the container, on the committed tree:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,248 passed, 0 failed, 14 ignored**, in 45.8 s, against 4,236 at
  `9aff9c13`.

**Eight mutations, eight killed**, two of which survived before the fixes:

| mutation | tests that fail | before the fixes |
|---|---|---|
| every numeric key dropped, the bases alone | 9 | **passed 7 of 8** |
| the second and third keys swapped | 2, the oracle among them | **the oracle passed** |
| the within-sample share key dropped | 6 | 3 |
| the samples-clearing key dropped | 3 | 1 |
| the cohort-read-total key dropped | 2 | 1 |
| the bases tie-break dropped | 2 | 2 |
| the share key ranks the smallest first | 7 | 4 |
| the bases tie-break reversed | 2 | 2 |

## 6. One thing worth keeping from how this review ran

**Every finding that mattered came from running something, and none from reading.** The oracle that
did not test its own claim, the bases silently agreeing with every numeric key, the invisible
argument transposition, and the 40 no-called samples were each produced by changing the code or the
fixture and watching what happened. The four agents agreed the comparator was correct — and it was.
What was wrong was what the tests could not see, and what the prose asserted about them.

**The prose was the weakest part of the step, not the code.** Four of the four numeric claims the
naming agent checked in the new doc comments were wrong, while every figure quoted *from* the design
documents was right. That is the same split the implementation skill records from an earlier plan in
this repo, and it is now two for two: claims about one's own fixture are the ones to re-run.
