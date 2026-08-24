# ng candidate alleles — C2: review and the fixes applied

*2026-08-24. Step C2 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md).
Reviewed at `1207018c` plus the step's working-tree diff, three agents in three isolated worktrees.*

---

## 1. Which categories ran

`reliability`, `naming`, and **design fidelity** — the third asked whether C3 could be built on this
step and told to answer by writing it.

## 2. What was actually wrong

### 2.1 No test could tell the specified ranking from production's — Blocker

Replacing the cap's comparator with production's key — the cohort read total, bases-tie-broken —
**left all 77 tests green.** The cause was in the fixtures rather than the assertions: every cap
test but one built its locus with one sample, and **at one sample the first ranking key and the
third are the same number divided by a constant**, so they induce the identical order. The one
multi-sample fixture gave all 400 alternatives the same share, the same clearing-sample count and
the same cohort total — a tie on every numeric key — and asserted only counts.

So the key spec §4.1 spends four paragraphs defending, and calls "not production's ranking, and the
difference is what it does at scale", was never the term that decided anywhere. Its worked case —
a real allele one sample carries at 30× scoring 15 reads against a mismapping artefact at 1 in 100
in 800 samples scoring 240 — is exactly the input the suite had none of. Getting it wrong cuts the
real allele, which step C3 turns into a missing genotype for its only carrier.

**Fixed** by the two tests the reviewer supplied, one for each regime spec §4.1 names: a
heterozygote at half of one sample's reads against a ten-sample artefact with more cohort reads,
and a three-reads-a-position locus where the shares tie and the sample count decides against a
larger cohort total.

### 2.2 The summary-to-bases pairing was never exercised — Blocker

Pairing each alternative's summary with its **neighbour's** bases also left all 77 green — the
exact mis-pairing `RankedAlternative` was invented at B2 to prevent, and whose own doc comment names
step C2 as the call site most at risk because every argument there is an index expression. Reversing
the bases comparison inside `compare_best_first` was equally invisible to `select_generic`.

The bases decide only at a numeric tie, and no cap test asserted **which** alleles survive a tie.
**Fixed** by a tie fixture where all three numeric keys tie and the bases alone separate the
alternatives — which is also the regime the cap actually meets at scale.

### 2.3 The one test at the cohort size the cap exists for asserted nothing the ranking could change — Major

The 400-private-allele locus asserted the verdict, the table length and the leftover length. All
three hold under **no ranking at all**. **Fixed** by asserting the identity of the survivors, which
at a total tie is a statement about the tie-break — and that turned the mis-pairing mutation from a
survivor into a failure.

### 2.4 Two fixtures listed their alternatives in descending read order — Major

So keeping the merge table's leading prefix satisfied them, and one asserted "the better-evidenced
of the two alternatives" in a message its fixture could not check — worse than saying nothing, since
a later reader will trust it. **Fixed** by reordering both so the best-evidenced alternative is not
the table's first, and asserting which survives.

### 2.5 No test reused a scratch across a truncated locus — Major

The cap sorts and truncates the buffer of surviving indices **in place**, so after a binding cap it
holds a short, reordered list that the next locus inherits. What makes that safe is `reset_for`'s
`clear`, which lives in the parent module and had nothing tying it to the in-place truncation this
step introduced. A stale index is either an out-of-range panic or an in-range index naming a
*different* allele of the new locus. **Fixed.**

### 2.6 Smaller things, all applied

- **`compare_best_first`'s `#[allow(dead_code)]` was stale** — C2 is the shipping caller its reason
  names. `allow` never warns when redundant, so nothing would have said so. Removed.
- `allowed` renamed `allowed_alternatives`: it is exactly the alleles-versus-alternatives
  distinction whose mutation killed six tests, and neither use site said which unit it counted.
- `selection_capped_at` now takes a `MaxCandidateAlleles` rather than a bare `u16`, so no call site
  has to respell in prose whether the number counts the reference.
- A test named `..._keeps_a_prefix_...` asserted a subset, and its own doc comment said subset.
- "the existing repeat-tract caller" was ng's own sibling in this file's vocabulary; it means
  production's, now said with the file.
- The `# Panics` paragraph claimed nothing else could panic on a wide locus, two lines from a
  `u32::try_from` that can. Reworded to name the width it refuses at.
- A cap of 12 in an earlier test was silently equal to the fixture's own 12-allele width, so it
  missed binding by zero; now 13.

## 3. Checked and found sound

- **C2 is spec §4.1 exactly on all four counts** — it truncates and never refuses, the reference
  never enters the ranked buffer, the ranking is `compare_best_first`, and `dropped` cannot include
  a rule-dropped allele because the buffer only ever held rule-clearing alternatives. The
  design-fidelity agent could not build an input where code and §4.1 disagree.
- **C3 builds on C2 unchanged**, and one of its findings simplifies C3: the plan says C3 must ask
  the rule again "over the alleles the cap cut", and it does not need to know which those were,
  because *this sample cleared the rule for a* implies *a was a candidate for the cap*. So the
  second count is `remap.candidate_for(allele).is_none() && reached_by(…)` — no cut list. **That
  also corrects my own dispatch note**, which said C3's count keys on the verdict: the *distinction*
  is load-bearing, the verdict field is a report.
- **The cap really does make a wide locus safe**: the reviewer ran a 70,000-allele locus at a cap of
  `u16::MAX` and got `Truncated { dropped: 4,465 }` with 65,535 candidates, one short of `admit`'s
  refusal.
- **Every number in the new prose checked out**, re-derived from the fixtures rather than read —
  the first time on this plan that has been true. The document figures hold at their sources too:
  HipSTR's `MAX_TOTAL_HAPLOTYPES = 1000`, production's repeat-tract `max_candidate_alleles: 24`,
  spec §4.1's "62 samples", spec §4.2's 23 of 53,935.

## 4. Raised, not applied — for Checkpoint C

1. **Spec §7 walks the bar across both ends of the committed range and never says what the *cap*
   does at one sample**, where truncating and refusing reach the same place. Measured by the
   reviewer: at one sample with ten alternatives at 3 reads each, nine clear the rule, five survive
   and the sole sample comes out missing.
2. The two items C1 already raised stand: arch §3.1's sentence about the order of the passes, and
   `CandidateAlleles::admit` accepting an empty allele where `new` refuses one.

## 5. Validation after the fixes

- `cargo fmt --check` clean; `cargo clippy --lib --tests --all-features -- -D warnings` clean;
  `cargo clippy --lib --all-features -- -D warnings` clean;
- `cargo test --lib` **4,276 passed, 0 failed, 14 ignored** in 37.6 s, against 4,265 at `1207018c`.

**Eight mutations, eight killed**, two of them survivors before the fixes:

| mutation | tests that fail | before |
|---|---|---|
| the cap ranks by production's cohort read total | 2 | **survived** |
| each summary paired with its neighbour's bases | 2 | **survived** |
| the ranking sort deleted | 5 | 3 |
| the kept prefix not put back into merge-table order | 1 | 1 |
| the cap boundary strict (`<` for `<=`) | 2 | 2 |
| the cap counting alleles where it counts alternatives | 6 | 6 |
| the cut count including rule-dropped alleles | 1 | 1 |
| the ranking reversed | 4 | 4 |

## 6. One thing worth keeping

**The single-sample fixture is this plan's recurring blind spot, and it has now cost four
Blockers.** B1's fixtures were shallow enough that the rule's floor decided and its share never
did; C1's were single-sample enough that a cohort sum matched the per-sample rule; C2's were
single-sample enough that a cohort total matched a within-sample share. Each time the term under
test was real, correct, and doing nothing that any fixture could observe.

The check that would have caught all four is one question asked of every fixture: **at this size,
is the term I am testing the one that decides — or would the simplest wrong rule give the same
answer here?**
