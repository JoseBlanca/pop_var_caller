# ng candidate alleles — C2: the cap

*2026-08-24. Step C2 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md),
branch `ng-candidate-alleles`, on top of `1207018c`. Design authority:
[`../../ng/spec/candidate_alleles.md`](../../ng/spec/candidate_alleles.md) §4, §4.1, §4.2 and
[`../../ng/arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md) §2.2, §2.5, §3.1.*

---

## 1. Plan

A locus is called over at most six sequences counting the reference. **Above that the list is cut
to the best and the locus is still called — never refused** (spec §4.1). Refusing is what HipSTR
does above 1,000 haplotypes and what production's repeat-tract path does above 24 candidates, and
at 63 accessions it costs 62 samples a locus they were called at perfectly well because one
accession carried something rare.

C2 inserts the cap between C1's two halves: rank the survivors with `compare_best_first`, keep what
the cap allows, and **sort the kept prefix back into merge-table index order**, because admission is
in table order and that is the order the VCF's `ALT` column comes out in.

## 2. Changes made

`src/ng/calling/allele_candidates/generic.rs`: the cap block inside `select_generic`, its doc
paragraphs, and eleven tests. `mod.rs`: `compare_best_first`'s now-stale `#[allow(dead_code)]`
removed — C2 is the shipping caller its reason named, and `allow` never warns when redundant.

**`Truncated { dropped }` counts only what the cap cut**, never what the admission rule dropped. A
reviewer building C3 corrected my understanding of why that matters: the *distinction* is
load-bearing, but C3 does not read the verdict to get it — a sample that cleared the rule for an
allele is by construction a sample that made it a candidate for the cap, so C3 asks
`remap.candidate_for(allele).is_none() && reached_by(…)` and needs no cut list.

## 3. What the review changed

Full account in [`../reviews/ng_candidate_alleles_c2_2026-08-24.md`](../reviews/ng_candidate_alleles_c2_2026-08-24.md).

**Two Blockers, both the same shape as the previous three steps': every cap fixture was one
sample.** At one sample the ranking's first key (the largest share of one sample's compared reads)
and its third (the cohort read total) are the same number over a constant — so **replacing the
ranking with production's cohort-total ranking left all 77 tests green.** That is the key spec §4.1
spends four paragraphs defending. And **pairing each summary with its neighbour's bases** was
equally invisible, because the bases decide only at a numeric tie and no fixture asserted which
alleles survive a tie — the exact mis-pairing `RankedAlternative` was built at B2 to prevent, at the
call site its own doc names as most at risk.

Three Majors: the 400-sample fixture asserted nothing the ranking could change; two fixtures listed
their alternatives in descending read order, so the merge table's own order gave the same answer;
and no test reused a scratch across a truncated locus, which is what the in-place sort-and-truncate
made new.

**Every number in the new prose was right this time** — the first step on this plan where the
naming reviewer's re-derivation from the fixtures found nothing.

## 4. Validation

- `cargo fmt --check` clean; both `clippy` gates clean;
- `cargo test --lib` **4,276 passed, 0 failed, 14 ignored** in 37.6 s, against 4,265 at `1207018c`.

**Eight mutations, eight killed**, two of them survivors before the review's fixes.

## 5. Follow-ups

- Spec §7 never says what the **cap** does at one sample, where truncating and refusing reach the
  same place. Raised at Checkpoint C.
- The two items C1 raised stand: arch §3.1's sentence about the order of the passes, and
  `CandidateAlleles::admit` accepting an empty allele where `new` refuses one.
