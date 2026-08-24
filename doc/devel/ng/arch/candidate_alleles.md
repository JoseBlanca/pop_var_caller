# ng step 6 — candidate alleles: types & interfaces

*Architecture draft, 2026-08-24. Companion to [`../spec/candidate_alleles.md`](../spec/candidate_alleles.md),
which argues every **why** below; this doc adds only the code shape. The repeat-tract path's
types are [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md), which extends these.*

*This document **is** step 6's real interface.*
[`ng_step_interfaces.md`](ng_step_interfaces.md) sketches the wire type in §2, the
`CandidateGenerator` recipe slot in §4, and marks the rename in §6 — all three deferring the
design; §5 below records what of that sketch survives. Shared
vocabulary — `CandidateAlleles`, `AlleleId`, `LocusKind` — is
[`calling_em_loop.md`](calling_em_loop.md) §2's and `calling/mod.rs`'s, and is used here, not
redefined.*

---

## Module home

`src/ng/calling/allele_candidates/`, the folder
[`module_layout.md`](module_layout.md) already reserves for step 6:

- `mod.rs` — everything both paths share: the config, the verdict, the leftover, the remapping,
  the output bundle, the ranking, and the fold that produces them.
- `generic.rs` — the ordinary path's entry point (§3.1).
- `ssr.rs` — the repeat-tract path's ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)).

**A folder rather than a file, and no trait.** Two paths with different inputs is two functions,
not two impls of one seam: they take different evidence, return different extras, and are chosen
by the locus kind rather than by a recipe. `ng_step_interfaces.md`'s `Box<dyn CandidateGenerator>`
does not survive — there is no bake-off here and nothing to swap.

---

## 1. What the module does, as a contract

**In:** one assembled `CohortObservation`
([`cohort_merge/build.rs:922`](../../../../src/ng/run/cohort_merge/build.rs)) and the run's
selection config. **Out:** a `LocusSelection` — the narrowed table, the verdict, the per-sample
leftover, and the remapping from the merge's allele indices to the new ids.

**Contract.**

- **Pure.** Reads one locus and the config; no other locus, no accumulated state, no clock, no
  random source. Same input, same bytes, at any worker count (spec §2, §8).
- **Infallible.** There is no error type. A locus that selects to the reference alone is a normal
  outcome (spec §6.2) and a locus above the cap is a normal outcome (spec §4.1). The remaining
  failure modes are caller bugs — a support row naming an allele outside the table, a non-finite
  `q_sum`, a sample with rows but no reads — and are assertions, the structural ones held in
  release, which is the convention [`read_likelihoods.md`](read_likelihoods.md) §1.1 sets for this
  module.
- **Allocation-free per locus** beyond the surviving table itself: the fold's buffers live in
  `CallingScratch` (§2.4).
- **The reference is never a candidate for removal** and stays at `AlleleId::REFERENCE`;
  `CandidateAlleles` makes that structural already
  ([`calling/mod.rs:101-118`](../../../../src/ng/calling/mod.rs)).

---

## 2. Types

### 2.1 The rule, and the constants

**The support rule is the merge's type, reused rather than copied**, so the two cannot drift
(spec §3). Only the default differs, and it gets its own named constant so a reader can see that
two rules share a shape and not a number.

```rust
/// The support one sample must lend one sequence for it to be called over, and the cap on
/// how many sequences a locus is called over at all.
pub struct CandidateSelectionConfig {
    /// `max(floor, ceil(share × that sample's reads at the locus))`, asked of each sample;
    /// one sample reaching it admits the sequence (spec §3).
    pub min_allele_support: MinAltReads,
    /// Counting the reference. Above it the list is cut to the best, never refused (spec §4).
    pub max_candidate_alleles: MaxCandidateAlleles,
}

/// How many alleles a locus may be called over, counting the reference. **At least two**:
/// below that the reference is the only survivor and every alternative becomes a
/// truncation, which is refusal under another name and is what spec §4.1 rules out.
/// `new` returns `None` below two, `new_or_panic` is its `const` path.
pub struct MaxCandidateAlleles(u16);

/// Floor 2 reads, share 10 in 100. **The floor is the merge's own and is defended there**
/// (`MinAltObs::DEFAULT`); the share is 10 in 100 where the merge's keep rule uses 2 (owner's
/// decision, 2026-08-24), set against a recall measurement rather than by it and unchanged
/// below 21 compared reads a sample (spec §3.3, §11 Q3).
pub const DEFAULT_MIN_ALLELE_SUPPORT: MinAltReads = MinAltReads {
    floor: MinAltObs::DEFAULT,
    share: MinAltReadShare::new_or_panic(0.10),
};

/// Six alleles including the reference — production's `DEFAULT_MAX_ALLELES_PER_RECORD`.
/// **Inherited, and measured to bind at about one tomato locus in 2,300 and none of the
/// human trio's** (spec §4.2). Soft.
pub const DEFAULT_MAX_CANDIDATE_ALLELES: MaxCandidateAlleles =
    MaxCandidateAlleles::new_or_panic(6);
```

**Three names here changed at Checkpoint A** (2026-08-24), on the code review's findings and
the owner's ruling. `support` became `min_allele_support` because in this crate "support"
already means the evidence a sample's reads showed (`AlleleSupportStats::num_obs`), and the
field holds a *threshold*; the constant moved with it. `max_candidate_alleles` became a
validated newtype because a bare `u16` let a cap of 0 or 1 be built and compiled — a reviewer
did exactly that — and neither value is a cap. And `MinAltReadShare::new_const` became
`new_or_panic`, because `pub const fn` is already visible in the signature while the `assert!`
is not.

**And the GATK half of the cap's lineage was wrong, inherited from production.** GATK's
`--max-alternate-alleles` defaults to **6 alternates** — `DEFAULT_MAX_ALTERNATE_ALLELES = 6`,
documented "Maximum number of alternate alleles to genotype"
([`GenotypeCalculationArgumentCollection.java:29`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypeCalculationArgumentCollection.java))
— so GATK genotypes over seven alleles. Production's constant counts a record's whole allele
set (`enforce_max_alleles` compares `unified.alleles.len()` and protects the reference from
pruning, [`per_group_merger.rs:1434`](../../../../src/var_calling/per_group_merger.rs)) and its
own doc comment nonetheless claims to match GATK's. **ng's cap is the tighter of the two by one
allele** — 21 genotypes against 28 at diploid — and the GATK clause is dropped from the
constant's documentation rather than repeated.

**`MinAltReads::reached_by` is the predicate, unwrapped** — it already takes a numerator and a
denominator ([`cohort_merge/mod.rs:466`](../../../../src/ng/run/cohort_merge/mod.rs)). Its
argument is named `non_reference_reads` because the merge asks it of a sample's pooled
non-reference reads; selection passes one allele's reads and a discovery round passes a narrower
count still (spec §3.4). **No wrapper**: a second spelling of one rule is how two rules become
different rules. *Impl-time: widen that method's doc comment to say the numerator is the caller's.*

### 2.2 The verdict

```rust
/// What selection did at one locus, beyond the list itself.
///
/// **There is no depth variant**, deliberately: depth is the merge's keep rule, asked once
/// upstream per sample, and production's version of it here is a cohort sum measured refusing
/// 98.6% of repeat tracts at one sample sequenced to 5× (spec §6.2).
#[non_exhaustive]
pub enum SelectionVerdict {
    /// Everything that cleared the bar is in the list. **Including when the list is the
    /// reference alone** — more than one built locus in four, on both benchmarks (spec §6.2).
    Selected,
    /// The cap bound and `dropped` alternatives were cut, lowest-ranked first (§2.5).
    ///
    /// **A `u32`, and it was a `u16`** (owner's decision, 2026-08-24). A review of step B2
    /// built a locus of 70,001 alternatives that all cleared the bar; narrowing the 69,995
    /// the cap cuts into a `u16` panics, so the one input on which this step would refuse a
    /// locus is the one spec §4.1's whole argument says must be truncated instead. Nothing
    /// upstream bounds a locus's allele table at 65,536 — `CandidateAlleles::admit`'s
    /// refusal at that width guards the *candidate* table, which the cap holds at six.
    Truncated { dropped: u32 },
    /// Repeat tracts only; never minted by `generic.rs`
    /// ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §3).
    NotPeriodic,
}
```

### 2.3 The leftover, and the remapping

```rust
/// One sample's reads whose sequence selection dropped, and the error mass they carry —
/// `read_likelihoods.md` §2.1's `unmatched_q_sum`, with the count it was missing.
///
/// **The second count is not decoration: it is what makes truncation defensible.** The mass is
/// the same under every genotype and cancels, so without it a sample whose true allele was cut
/// is scored confidently against a set that does not contain it and an invented genotype comes
/// out (spec §4.1, §5). Drop it and refusing the locus becomes the correct policy again.
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub struct UnmatchedSupport {
    pub num_reads: u32,
    /// Σ `ln P(error)` over them — summed straight from the merge's own per-row `q_sum`,
    /// so it is never re-derived. Zero (not negative) where nothing was dropped.
    pub q_sum: f64,
    /// Of those reads, the ones on an allele that cleared the bar **for this sample** and was
    /// then cut by the cap. **Non-zero means this sample's genotype is emitted as missing.**
    pub earned_reads_cut_by_the_cap: u32,
}

impl UnmatchedSupport {
    /// `earned_reads_cut_by_the_cap > 0`, named so the rule is greppable.
    pub fn genotype_must_be_missing(&self) -> bool;
}
```

**Why the second count and not `num_reads > 0`** (owner's decision, 2026-08-24). The bar drops
alleles almost nobody showed — 13,166 of 15,474 alternatives on the GIAB trio at 300×, spec §3.3 —
so nearly every sample has a non-zero pool at nearly every locus and a rule keyed on it would
emit a missing genotype almost everywhere. The cap only ever cuts alleles that cleared the bar for
*somebody*; asking whether it cleared for *this* sample is what makes the rule fire exactly where
a real allele was taken from that sample.

```rust

/// For each allele of the merge's table, in that table's own index order: the id it now has
/// among the candidates, or nothing where selection dropped it.
///
/// **The evidence builder cannot be written without this.** `CandidateAlleles` ids are dense
/// and in admission order ([`calling/mod.rs:170-186`](../../../../src/ng/calling/mod.rs)),
/// while `SupportedAllele::allele` indexes the merge's table
/// ([`build.rs:1089`](../../../../src/ng/run/cohort_merge/build.rs)); after narrowing the two
/// are different numbers and nothing else records the correspondence.
pub struct AlleleRemap {
    to_candidate: Box<[Option<AlleleId>]>,
    /// How many have been admitted — the id the next admission must carry.
    num_admitted: u32,
}

impl AlleleRemap {
    /// Every allele dropped; selection fills the survivors in as it admits them.
    pub fn with_all_dropped(table_len: usize) -> Self;

    /// The candidate id for one of the merge table's alleles, or `None` where it was
    /// dropped. `None` is "absent", never a sentinel id. Asserts on an index outside the
    /// merge's table, which is a support row naming an allele that does not exist.
    pub fn candidate_for(&self, table_index: usize) -> Option<AlleleId>;

    /// Record that `table_index` survived as `candidate`. Asserts three ways.
    pub fn admit(&mut self, table_index: usize, candidate: AlleleId);

    pub fn table_len(&self) -> usize;
    pub fn num_admitted(&self) -> usize;
}
```

**The count is carried, and it is what closes the hole a bounds check cannot see.** An earlier
draft of this section declared `candidate_for` alone; nothing could then build a remapping, and
the writer that had to be added guarded the merge index three ways and the candidate id not at
all. Two reviewers independently compiled `admit(1, AlleleId(1))` followed by
`admit(2, AlleleId(1))`: both indices in range, each written once, and two of the merge's alleles
resolving to one candidate. The evidence hand-off of §3.2 re-keys the merge's support rows through
this map, so that state hands two different sequences' reads to one candidate and the read
likelihood scores two alleles as one — with an ordinary-looking genotype coming out.

`admit` therefore asserts that `candidate` is **the next dense id**, which is what
`CandidateAlleles::admit` returns by construction. **The check relates the id to the admission
count and never to `table_index`**, so it is independent of the order a caller walks the table in
— which matters because the repeat-tract path admits by ladder rung rather than in table order
([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §3.1). It also makes `num_admitted` a
field read rather than a scan.

**The merge-table index stays a bare `usize` across this surface**, considered and kept at
Checkpoint A. A newtype would make the two index spaces impossible to confuse, but they are
already different *types* — a candidate is an `AlleleId`, a merge index is a `usize`, and turning
one into the other takes a deliberate cast. The upstream producer `SupportedAllele::allele` is a
bare `usize` too, so a newtype here alone would add a wrap at every call site and protect against
a confusion the type system already refuses.

### 2.4 What one call returns, and where its buffers live

```rust
/// Everything selection produces at one locus.
///
/// **The fields are private and `new` is the only door.** It asserts both invariants below;
/// public fields would have made that check bypassable by a struct literal, on a value whose
/// defect is a wrong genotype rather than a crash. Read through `alleles()`, `verdict()`,
/// `unmatched()` and `remap()`, or take ownership with `into_parts()`.
pub struct LocusSelection {
    alleles: CandidateAlleles,
    verdict: SelectionVerdict,
    /// **Parallel to `CohortObservation::per_sample`** — same length, same order, so entry `i`
    /// belongs to that entry's sample. Not indexed by the run's sample order, because
    /// `per_sample` holds only the covering samples and each names its own
    /// ([`build.rs:929-941`](../../../../src/ng/run/cohort_merge/build.rs)).
    unmatched: Vec<UnmatchedSupport>,
    remap: AlleleRemap,
}

impl LocusSelection {
    /// Asserts that `unmatched` is one entry per covering sample, and that the remapping
    /// admitted exactly as many alleles as the table holds.
    pub fn new(
        alleles: CandidateAlleles,
        verdict: SelectionVerdict,
        unmatched: Vec<UnmatchedSupport>,
        remap: AlleleRemap,
        covering_samples: usize,
    ) -> Self;
}

/// The fold's buffers, one set per worker. **A nested section of `CallingScratch`**
/// ([`calling_em_loop.md`](calling_em_loop.md) §2) rather than a scratch of its own: the same
/// worker runs selection and then the loop on the same locus, so a second per-worker
/// allocation buys nothing.
///
/// **`reset_for` destructures rather than naming its fields through `self`**, so that a
/// buffer added later — the repeat-tract path commits to adding one
/// ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §5) — fails to compile there
/// instead of silently carrying one locus's values into the next.
pub struct SelectionScratch {
    per_allele: Vec<AlleleSummary>,          // cleared and refilled per locus
    ranked_table_indices: Vec<u32>,          // merge table indices, sorted by §2.5
}
```

`AlleleSummary` is private to the module — one allele's fold across samples: the largest share of
one sample's reads it took, how many samples cleared the bar, and its cohort read total.

**Two of its three fields are readable from outside as scalars**, through
`SelectionScratch::best_within_sample_share_of` and `SelectionScratch::cohort_reads_of` (step D1).
They exist for `examples/ng_candidate_selection_probe.rs`, which reports the ranking's own keys and
would otherwise have to recompute them — the duplicate rule that step deletes. **Two scalars rather
than the type**, so the shape of the computation stays inside the module and nothing outside can
hold or build a summary; the third field has no reader and is not exported.

**It does *not* carry the reads and mass it would contribute to the leftover**, which an earlier
draft of this section asked for. Three reasons, and the third is the one that settles it:
`AlleleSummary` is per allele with no sample axis while `LocusSelection::unmatched` is per
covering sample, so no per-allele total can produce the output; survival is not known during the
pass that fills the summary, so the number could not be filled anyway; and **a cohort total is a
sum in allele-major order where the leftover's own oracle demands the per-sample rows' sum**, so
the bitwise check that step would be isolated for would fail by construction if the total were its
source. Nothing reads it — the ranking uses the share, the samples clearing and the cohort read
total — and §6's future `q_sum` bar reserves nothing here either. *Corrected at Checkpoint A,
2026-08-24, on a finding two reviewers reached independently;
[`../impl_plan/candidate_alleles.md`](../impl_plan/candidate_alleles.md) step B2 carries the same
correction.*

### 2.5 The ranking

```rust
/// One alternative as the ranking reads it: the fold's summary and the bases that break the
/// last tie, travelling together.
struct RankedAlternative<'bases> {
    summary: AlleleSummary,
    bases: &'bases [u8],
}

/// Order two alternatives best-first for the cap: the largest share of one sample's reads
/// first, then how many samples cleared the bar, then the cohort's read total, then the bases
/// (spec §4.1 — why this and not production's cohort read total).
///
/// **The name carries the direction because the return cannot.** `Ordering::Less` means `left`
/// belongs earlier in a best-first list; `min_by` gives the best allele and **`max_by` gives
/// the worst**.
///
/// **The bases are the tie-break that cannot tie**, which is what makes the order independent
/// of the merge's own allele order rather than inheriting it. A three-way tie before that point
/// is rare enough that the byte comparison costs nothing.
///
/// **Shares compare with `f64::total_cmp`** — a total order, so there is no NaN branch and no
/// partial-order footgun. A `NaN` share is a caller bug and asserts in the fold, not here.
fn compare_best_first(left: RankedAlternative<'_>, right: RankedAlternative<'_>) -> Ordering;
```

*Corrected 2026-08-24, on step B2's review, and both halves were measured rather than argued.*
**The name was `ranks_above`**, which is the shape Rust reserves for `-> bool` predicates
answered with a value whose obvious reading is the opposite — and the cost is concrete:
`max_by(ranks_above)` returns the *worst*-ranked allele, compiles, and is what somebody reaching
for "the best one" would write. **The signature took four positional arguments** — two summaries
and two base slices — and a reviewer swapped the two slices at a call site: it compiled, `clippy`
was silent, and the test still passed, because the bases only decide when all three numeric keys
tie. **The mis-pairing is therefore invisible at exactly the loci where the ranking does its
work**, and the shipping caller is worse than any test: step C2 sorts a buffer of table indices,
so every argument is an index expression.

---

## 3. The interface

### 3.1 The ordinary path

```rust
/// Narrow one locus's allele table to the sequences worth calling over (spec §3, §4).
///
/// **One pass over the locus's rows.** For each sample: its reads at the locus are the sum of
/// its rows, which is that sample's compared reads because the merge admits only complete
/// observations onto alleles (spec §1.3); each row then feeds its allele's summary. A second
/// pass admits the survivors in table order, applies the cap, and fills the leftover.
///
/// The reference is admitted first and is exempt from both the bar and the cap.
pub fn select_generic(
    observation: &CohortObservation,
    config: &CandidateSelectionConfig,
    scratch: &mut SelectionScratch,
) -> LocusSelection;
```

**Read-group rows are pooled here, and that is the one place pooling is correct.** The merge keys
its rows on `(allele, read group)` because the likelihood must not pool them
([`build.rs:1089`](../../../../src/ng/run/cohort_merge/build.rs)); the bar counts reads, and a read
is a read whichever lane it came from, so it sums over groups —
`SampleSupport::pooled_support_for` is the existing method that does it
([`:1070`](../../../../src/ng/run/cohort_merge/build.rs)).

**Partial observations are not read.** They are held on their own axis with their own bases and
witnessed stretch ([`PartialObservation`, `build.rs:1130`](../../../../src/ng/run/cohort_merge/build.rs)),
they count toward no bar, and they do not enter the leftover (spec §5.1).

### 3.2 What the caller does next

Selection does not build `GenericSampleEvidence`; the loop's input edge does
([`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md)). What this module owes it is
`remap` plus `unmatched`, and the shape of the hand-off is fixed:

```text
for each covering sample i:
    GenericSampleEvidence {
        supported:       rows of per_sample[i] whose remap.candidate_for(row.allele) is Some,
                         re-keyed to that AlleleId,
        unmatched_q_sum: unmatched[i].q_sum,
        partials:        per_sample[i].partials, untouched,
    }
```

**`GenericSampleEvidence` gains a field.** `read_likelihoods.md` §2.1 declares `unmatched_q_sum`
alone; the count of §2.3 has to travel with it or nothing downstream can tell a truncated sample
from a clean one. *That is an edit to that document and is not made here.*

---

## 4. Design decisions — decided

- **No trait, two functions.** Two paths with different inputs and different extras, chosen by
  the locus kind; `ng_step_interfaces.md`'s `Box<dyn CandidateGenerator>` had no second impl to
  swap in. — spec §1, module home above.
- **The support rule is `MinAltReads`, reused.** One home for the contract so the merge's rule and
  the allele rule cannot become different rules. — spec §3.
- **The share is 10 in 100 and the floor is 2.** The floor is the expensive knob: at 30×, raising
  it to 3 loses five true alleles to keep 1,539 alternatives where the share at 10 in 100 loses
  two to keep 1,601. The share's own two are the price the owner set it at, against the recall
  measurement rather than by it. — spec §3.3, §11 Q3.
- **Truncate, ranked by the largest within-sample share.** Production ranks by the cohort read
  total and truncates private alleles first at scale; refusal loses the good alleles with the bad.
  — spec §4.1.
- **The leftover carries a read count.** Truncation and the count stand or fall together. — spec §5.
- **The verdict has no depth variant.** — spec §6.2.
- **Selection returns a remapping.** Dense candidate ids and merge table indices are different
  numbers after narrowing, and nothing else records the correspondence. — this doc §2.3; the spec
  does not name it because it is a code shape, not a design choice.
- **`SelectionScratch` nests inside `CallingScratch`.** One allocation per worker, not two. —
  [`calling_em_loop.md`](calling_em_loop.md) §2's convention.
- **Infallible, assertions for caller bugs.** — spec §8, matching
  [`read_likelihoods.md`](read_likelihoods.md) §1.1.
- **`DEFAULT_MAX_CANDIDATE_ALLELES` lives here and `CallingLoopConfig` reads it.**
  [`calling_em_loop.md`](calling_em_loop.md) §2.1 declares the same constant because the loop
  enforces the cap after a discovery round; **it must be the same constant, not a second one with
  the same value.** — spec §6.3.

---

## 5. Reconciliation with existing code

| this doc's name | existing code | how it converges |
|---|---|---|
| the support rule | [`MinAltReads`, `cohort_merge/mod.rs:424`](../../../../src/ng/run/cohort_merge/mod.rs); `required_of` [`:451`](../../../../src/ng/run/cohort_merge/mod.rs), `reached_by` [`:466`](../../../../src/ng/run/cohort_merge/mod.rs) | **the type itself, unchanged.** Only `DEFAULT_MIN_ALLELE_SUPPORT` is new, and it reuses `MinAltObs::DEFAULT` ([`:317`](../../../../src/ng/run/cohort_merge/mod.rs)) |
| the input | [`CohortObservation`, `build.rs:922`](../../../../src/ng/run/cohort_merge/build.rs); `SampleSupport` [`:965`](../../../../src/ng/run/cohort_merge/build.rs), `SupportedAllele` [`:1089`](../../../../src/ng/run/cohort_merge/build.rs), `AlleleSupport` [`:1245`](../../../../src/ng/run/cohort_merge/build.rs) | read, unchanged. `q_sum` on `AlleleSupport` is the leftover's only source |
| pooling a sample's rows for one allele | [`SampleSupport::pooled_support_for`, `:1070`](../../../../src/ng/run/cohort_merge/build.rs) | called as-is (§3.1) |
| the output table | [`CandidateAlleles`, `calling/mod.rs:86`](../../../../src/ng/calling/mod.rs); `new` [`:101`](../../../../src/ng/calling/mod.rs), `admit` [`:184`](../../../../src/ng/calling/mod.rs) | built and merged; selection seeds with `new` and fills with `admit`. **No new allele-table type is minted** |
| the id | [`AlleleId`, `types.rs:304`](../../../../src/ng/types.rs), `AlleleId::REFERENCE` | reused; `AlleleRemap` holds `Option<AlleleId>` and never a sentinel |
| the cap's value | [`DEFAULT_MAX_ALLELES_PER_RECORD`, `per_group_merger.rs:57`](../../../../src/var_calling/per_group_merger.rs) | the number is inherited and declared inherited; **the ranking is not ported** — production's is [`enforce_max_alleles`, `:1434`](../../../../src/var_calling/per_group_merger.rs) |
| the leftover | [`pool_dropped_other_scalars`, `per_group_merger.rs:1555`](../../../../src/var_calling/per_group_merger.rs), consumed [`:1979`](../../../../src/var_calling/per_group_merger.rs) | same shape, plus the read count. **Do not repeat its contamination-path omission** ([`posterior_engine.rs:2429`](../../../../src/var_calling/posterior_engine.rs)) |
| the consumer's field | `GenericSampleEvidence::unmatched_q_sum`, [`read_likelihoods.md`](read_likelihoods.md) §2.1 | filled by `unmatched[i].q_sum`; that struct gains the read count (§3.2) |
| the superseded sketch | `AlleleCandidates` + `Admission` ([`ng_step_interfaces.md`](ng_step_interfaces.md) §2), `CandidateGenerator` (§4), the rename note (§6) | **superseded.** `AlleleCandidates` was already renamed to `CandidateAlleles`; `Admission`'s four variants become `SelectionVerdict`'s three (§2.2); the trait is dropped |

---

## 6. Open items

- `OPEN:` **whether a `q_sum` bar joins the read-count bar at high depth** — spec Q1's neighbour,
  recorded in spec §3.3 with what would trigger it. Nothing here reserves a field for it: the
  merge already carries `q_sum`, so adding it later is a config field and a term in the fold.
- `OPEN:` **whether the cap becomes load-bearing past 63 samples, and whether the ranking then
  matters** — spec Q2. No code shape depends on the answer; both are already switchable values.
- *Impl-time, done:* `MinAltReadShare` gained `new_or_panic`, a `const` path beside its fallible
  `new`, because `DEFAULT_MIN_ALLELE_SUPPORT` cannot reach the type's private field. Both
  constructors now share one `const fn is_a_fraction_of_one`, so they cannot come to disagree
  about what a legal share is — they had two copies of the test, and a review's mutation of one
  left the other's tests green.
- *Impl-time:* whether `AlleleRemap` is `Box<[Option<AlleleId>]>` or a bitset plus a prefix sum.
  The first is written here because a locus holds a handful of alleles; measure before changing.
  **What a locus allocates, for that measurement's baseline:** the surviving table, one
  `Box<[Option<AlleleId>]>` as long as the *merge's* table, and one `Vec<UnmatchedSupport>` one
  entry per covering sample. §1's "allocation-free per locus" is about the fold's working
  buffers, which live in the scratch; these three are the output.

---

## Test & bench shape

Tests live in `#[cfg(test)] mod tests` beside each file, per the repo rule. Spec §12 lists what
they assert; the two that pin *this* document's shapes rather than the spec's rules:

- **the remapping** — after a narrowing that drops the middle allele of five, every surviving
  row's `candidate_for` returns the dense id and the dropped one returns `None`; feeding the
  result through §3.2's hand-off reproduces the evidence view exactly.
- **the leftover is the merge's own arithmetic** — the pool equals the sum of the dropped rows'
  `q_sum` to the last bit, not a re-derivation from counts and a rate.

The regression anchor is spec §12's last entry: the GIAB trio and the tomato panel through the
calling loop with real candidates, which is the blocker
[`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md) records. The measurement harness
already exists — `examples/ng_candidate_selection_probe.rs` — and **calls this module rather than
carrying a copy of it** (step D1, 2026-08-24), so the numbers the spec quotes are the shipped
code's. It reproduced them: 4,177 and 7,478 built loci on the trio at 30× and 300×, 53,935 on the
tomato panel, every bar's kept-alternative total, every cap's binding count, and the leftover at
about 4 tomato reads in every 1,000 (143,712 of 39,589,086). **One figure changed with the code** —
the two cap rankings keep different alleles at 19 tomato loci where the standalone copy found 17,
because the within-sample share is now maximised over the samples that cleared the bar. **One was a
slip in a table**, where the widest locus at 16 samples asked carries 10 alternatives and was
written as 14. Spec §4.2 records both.
