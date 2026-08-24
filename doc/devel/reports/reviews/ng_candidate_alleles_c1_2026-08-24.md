# ng candidate alleles — C1: review and the fixes applied

*2026-08-24. Step C1 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md).
Reviewed at `eff6cf16` plus the step's working-tree diff, three agents in three isolated worktrees.
Implementation report: [`../implementations/ng_candidate_alleles_c1_2026-08-24.md`](../implementations/ng_candidate_alleles_c1_2026-08-24.md).*

---

## 1. Which categories ran

`reliability`, `naming`, and **design fidelity** — the third asked whether C2 and C3 could be built
on this step, and told to answer by writing them. Three rather than the six of B1: the step is one
public function of forty lines with no fallible path, no concurrency and no new public type, and
the three categories that have produced every finding on this plan so far are these.

## 2. What was actually wrong

### 2.1 The admission rule could be replaced by a cohort read total — Blocker

Spec §3.2 is the rule the whole module hangs on: **no term of the admission bar may read the
cohort**, because otherwise a sample's candidate list depends on who else is in the run. Replacing
`cleared_the_bar()` with `cohort_reads >= 2` **left all 65 tests green.**

The two rules separate only when two or more samples each lend an alternative *less* than the
floor and their reads pool over it, and no fixture built that: every fixture showing an alternative
showed it from one sample. The nearest one adds a second sample showing **reference reads only**,
which adds nothing to any alternative's cohort total — so it moves a cohort *denominator* and not a
cohort *numerator*. **Its own doc comment claimed it was "the test that fails first if a cohort
term ever creeps into the admission rule."** Measured, it does that for a
majority-of-covering-samples rule and not for a pooled count.

**Fixed** by the fixture the file had no version of: two samples lending one read each, pooling to
the floor of two, neither reaching it alone. And the false claim on its neighbour corrected to say
which half it covers.

### 2.2 The rule's share could be dropped — Blocker

Keeping the configured floor and forcing the share to zero also left all 65 green — and unlike
passing the shipped default, it keeps `config` used, so `-D warnings` sees nothing either. Every
fixture was built at a share of 0.0 except one, where `ceil(0.5 × 3) = 2` is the floor.

**The branch no test exercised is the one the high-depth benchmark spends all its time in.** The
shipped share of 5 in 100 binds above 41 compared reads; the GIAB trio runs at 30× and 300×. A rule
silently degraded to its floor there admits sequencing error as a candidate — 10 reads in 300 is
about the error rate — and the symptom is a longer `ALT` list, not a crash.

**Fixed** by a deep fixture: 10 alternative reads of 300 compared, where `ceil(0.05 × 300) = 15`
refuses what a floor of 2 would admit.

### 2.3 A covering sample with no support rows was in no fixture — Major

`covering_samples` sizes the leftover *and* is the length `LocusSelection::new` checks against, so
a wrong count is self-consistent and the constructor cannot see it. Counting only the samples that
have rows survived all 65 tests — the three-sample fixture's samples all have rows.

The merge builds the case it misses: a sample that covered the locus and whose reads all stopped
inside it has partials and nothing else. **The cost is not a short vector but a shifted one** — the
leftover is parallel to `per_sample` by position, so every later sample's leftover slides onto its
neighbour, and once C3 fills these that is a missing genotype for a sample that lost nothing and an
invented one for the sample that did.

**Fixed**; the fixture helper already existed in the shared module and simply had no caller here.

### 2.4 Nothing asserted the leftover was zeroed — Major

Every test read `.len()` or `.is_empty()`. Filling the vector with a non-default
`UnmatchedSupport` left all 65 green — and `earned_reads_cut_by_the_cap > 0` is
`genotype_must_be_missing()`, so a non-zero third field on that line **emits every covering sample
at every locus as a missing genotype, the whole run silently no-called.** That line is precisely
what step C3 rewrites, which is when the slip becomes likely.

**Fixed** by asserting the value beside the length.

### 2.5 A merge table holding one sequence twice was admitted twice — Major

Probed on the clean tree: both copies become candidates and both merge indices remap onto distinct
ids. Two identical `ALT` sequences reach the VCF and the read likelihood scores one sequence as two
alleles, splitting its evidence. `CohortObservation::alleles` documents distinctness and
`AlleleTable` enforces it by interning on the bytes, so this is a producer invariant rather than an
input — **held now by a `debug_assert!`**, not a release assertion, because the check is a scan of
the table.

### 2.6 A wrong mechanism in the step's own prose — Minor, and the one worth naming

The `# Panics` note said the empty-table assertion's alternative was `CandidateAlleles::new`
panicking about an empty reference allele. The reviewer deleted the assertion and ran it:
`index out of bounds: the len is 0 but the index is 0`, two statements earlier, where the
reference's bases are read. `CandidateAlleles::new` is never entered, and its own assertion is on
empty *bases*, which an empty table never produces.

**A wrong mechanism is worse than a wrong number**, because it sends the next reader hunting a
symptom that does not occur. Corrected with what actually happens, and marked as measured.

### 2.7 Smaller things, all applied

- The rule was written `max(2 reads, share × compared reads)` where the spec and `required_of` both
  have `ceil(share × …)`.
- The paragraph explaining a `#[allow(dead_code)]` survived the attribute's removal, sending the
  reader looking for something no longer there.
- `LocusSelection::new`'s doc still said "the fields stay public because arch §2.4 declares them
  so", contradicting the type's own doc two paragraphs above and the code since Checkpoint A.
- The lifted fixtures spelled `pub(in crate::ng::calling::allele_candidates)` seven times where
  `pub(super)` means the same thing.
- A test for the scratch reused **narrow then wider**, the direction its sibling does not walk.

## 3. Checked and found sound

- **Fidelity to spec §3, §6.1 and §6.2 is exact.** The design-fidelity agent re-ran one of the
  author's own mutations rather than taking it on trust and got the three failures reported.
- **C2 and C3 build on C1 without a single change to it** — not a signature, not a buffer, not a
  type. The agent wrote both in its worktree: C2 is 27 lines inside `select_generic`, C3 a 30-line
  free function plus a three-line `map`. Library suite with both added: 4,269 passed.
- **C3's hardest question is reachable, and more cheaply than the plan expects.** The leftover needs
  *neither* buffer of the scratch — only the sample's rows, the finished remapping and the config.
  And the cap-versus-bar distinction needs no test at all: a sample that cleared the rule for an
  allele is by construction a sample that put it into the survivor list, so an allele *this* sample
  earned and the remapping does not hold **can only have been cut by the cap**.
- **The `ranked_table_indices` contract is safe with both of C2's sorts in place.** `reset_for`
  clears and reserves unconditionally, and sorting and truncating never touch capacity.
- **The round-trip test's hole in the middle is load-bearing, and `AlleleRemap`'s own assertions are
  not doing its work** — a counting remap satisfies all three of them and is caught by the test's
  own assertion.
- Every design-document figure quoted in the new prose checked out at source, and the reviewer
  re-derived all eight fixtures' arithmetic from the fixtures themselves.

## 4. Raised, not applied — for Checkpoint C

1. **The cap's by-catch is samples, not alleles, and at scale it is most of them.** At the case the
   cap exists for — 400 samples each carrying a different private allele — the verdict is
   `Truncated { dropped: 395 }` and **395 of the 400 samples are emitted as missing**, because each
   earned the allele the cap took away. That is the design behaving exactly as specified (spec
   §4.1, §5's second count). But spec §4.1's reassurance that "only the samples that earned a cut
   allele are affected" is measured at 63 tomato accessions where the cap binds at 23 of 53,935
   loci; **at several hundred samples the same sentence means almost everybody.** This is spec §11's
   Q2 made concrete rather than a defect in any step.
2. **Arch §3.1's sentence about the order of the passes cannot be implemented as written** — the
   cap must precede admission and the leftover must follow it. Recorded in the code's own doc
   comment; the document is the owner's to edit.
3. **`CandidateAlleles::admit` accepts an empty allele where `CandidateAlleles::new` refuses one.**
   `new` refuses it because an empty `REF` would reach the VCF as an unparseable record; `admit`
   has the identical hazard one column to the right and no check. The guard belongs in
   `calling/mod.rs`, which this plan does not own.
4. **C2 must carry a test pinning admission order at a binding cap.** The agent deleted the
   sort-back that returns the kept prefix to merge-table order and no C1 test noticed, though the
   whole `ALT` column was permuted.

## 5. Validation after the fixes

All in the container, on the committed tree:

- `cargo fmt --check` clean; `cargo clippy --lib --tests --all-features -- -D warnings` clean;
  `cargo clippy --lib --all-features -- -D warnings` clean;
- `cargo test --lib` **4,265 passed, 0 failed, 14 ignored** in 45.6 s, against 4,249 at `eff6cf16`.

**Nine mutations, nine killed**, four of them survivors before the fixes:

| mutation | tests that fail | before |
|---|---|---|
| admission keyed on the cohort read total | 2 | **survived** |
| the configured share dropped, the floor kept | 1 | **survived** |
| only samples with support rows counted as covering | 1 | **survived** |
| the leftover filled non-zero | 1 | **survived** |
| the admission rule inverted | 6 | 6 |
| the remapping off by one | 7 | 7 |
| the reference offered to the rule | 8 | 8 |
| the survivors admitted in reverse table order | 3 | 3 |
| every survivor admitted with the reference's bases | 3 | 3 |

## 6. One thing worth keeping from how this review ran

**Three steps, three times two Blockers, and all six were tests that could not fail.** Not one of
them was wrong code — the fold, the ranking and the admission pass each computed the right answer
on every input any agent could build. What was wrong each time was that a property nobody had
questioned was pinned by nothing, and the only way to see it was to break the property and watch
the suite stay green.

The specific shape recurs too: **a fixture that exercises a rule at a depth or a cohort size where
the rule's two halves agree.** B1's fixtures were all shallow enough that the floor decided; C1's
were all shallow enough that the floor decided *and* single-sample enough that a cohort sum
matched. The lesson is cheap to state and apparently expensive to remember: **a fixture must be
built at the size where the term under test is the one that decides**, and raising a threshold in
the fixture is not the same as raising the depth.
