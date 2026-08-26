# ng step 6 at a repeat tract — types & interfaces


> **⚠ 2026-08-26 — the repeat tract's genotype prior no longer indexes its seed by the
> cohort's modal repeat count**, so every sentence below that justifies carrying the mode
> *for the prior* has lost its reason. The seed is now the fitted **length spectrum** the
> joint repeat fit produces per stratum, indexed by whole-repeat offset from the
> **reference** tract length, which every locus already knows
> ([`../spec/population_diversity.md`](../spec/population_diversity.md) §4.2; built by step
> E2e of [`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md)). Handing the mode
> where the reference length belongs is now a measured error: it moves 0.595 of the prior's
> mass off the reference length onto 0.091 on `seed_ssr`'s own fixture. **Whether selection
> should still carry the mode is undecided and is this document's to settle** — it may have
> uses of its own; what it no longer has is this one.
*Architecture draft, 2026-08-24. Companion to
[`../spec/candidate_alleles_ssr.md`](../spec/candidate_alleles_ssr.md), which argues every **why**
below. **Extends [`candidate_alleles.md`](candidate_alleles.md)** — the config, the verdict, the
leftover, the remapping, the output bundle and the ranking are that document's and are used here,
not redefined. Read it first.*

*Blocked on one field the merge does not carry (§1). The types below are settled; the entry point
has no producer for one of its arguments until that lands.*

---

## Module home

`src/ng/calling/allele_candidates/ssr.rs`, beside `generic.rs` under the shared `mod.rs`
([`candidate_alleles.md`](candidate_alleles.md), *Module home*). One file, no trait: this is the
second of the two path functions, not an impl of a seam.

---

## 1. The blocker, as a signature

**`CohortObservation` carries no locus kind**, so nothing in it names the motif, and without a
motif there is no period and no repeat count
([`build.rs:922-936`](../../../../src/ng/run/cohort_merge/build.rs); spec §2). The entry point
therefore takes the detail as a separate argument, which is the honest shape today and the shape
that survives when the merge fixes it:

```rust
// OPEN (blocked): `detail` has no producer. `SampleLocusObservations` carries
// `LocusKind::Ssr(SsrDetail)` (`locus_generation/mod.rs`), the merge's closer branches on that
// kind (`cohort_merge/close.rs:110-125`), and the assembler drops it. The merge owes
// `pub kind: LocusKind` on `CohortObservation`; when it lands, `detail` comes off it and this
// argument goes away.
```

Nothing else in this document waits on it.

---

## 2. Types

### 2.1 The ladder

**A rung is a repeat count and the sequences observed at it.** The ladder is built once per locus
from the merge's table plus the motif, and it is what nomination reads.

```rust
/// One locus's observed tract sequences, grouped by how many repeats they carry.
///
/// The key is `bases.len() / motif.period()`, floored — **the same integer the genotype prior
/// indexes its geometric seed by** (`calling_priors.md` §5), which is the reason the grouping
/// exists at all and the reason both must compute it one way. Sequences inside a rung stay
/// distinct: the grouping is for nomination and for the prior, never a merge of evidence
/// (spec §3).
pub struct RepeatLadder {
    /// Ascending by repeat count. Each entry holds the merge-table indices of the sequences
    /// at that count, ascending, so a rung never owns bases.
    rungs: Vec<Rung>,
    /// The cohort's modal repeat count — the most-supported rung, ties to the shorter. The
    /// prior's decay centre (`calling_priors.md` §5), and **not** the reference's count.
    modal_repeat_count: u32,
}

struct Rung { repeat_count: u32, table_indices: Vec<u32>, cohort_reads: u64 }
```

`RepeatLadder` is public only so far as the tests need it; the entry point returns what the
consumers want (§2.3) rather than the ladder itself.

### 2.2 Configuration

```rust
/// The repeat-tract path's own settings, beside the shared ones.
pub struct SsrSelectionConfig {
    /// The support rule and the cap — `candidate_alleles.md` §2.1's type, reused. The share
    /// is 5 in 100 here too, which is why it is not a second field (spec §5).
    ///
    /// **The cap is not the sibling's.** This path builds it from
    /// [`DEFAULT_MAX_CANDIDATE_ALLELES_SSR`] — **32**, against the ordinary path's six —
    /// because a tract genuinely carries more alleles than a SNP does and HipSTR, the
    /// nearest comparator, caps the count at nothing at all (spec §12, Q2, settled
    /// 2026-08-24). The type is shared; only the number differs, and it costs 528 diploid
    /// genotypes a sample a locus against six's 21.
    pub shared: CandidateSelectionConfig,
    /// The share of one sample's spanning reads that may sit off the motif grid before that
    /// sample is judged non-periodic (§3).
    pub max_off_grid_share: f64,
    /// Copies per genome — how many rungs a sample promotes (§3.1). From
    /// `FrozenParameters::ploidy`, never a constant here.
    pub ploidy: Ploidy,
}

/// One read in ten off the motif grid. **Production's `max_out_of_frame_frac`, inherited and
/// never measured — by them or by us** (`candidate_set.rs:85`; spec Q1). Soft: the measurements
/// in spec §4.1 and §5 were all taken with this gate switched off, so nothing there constrains
/// it.
pub const DEFAULT_MAX_OFF_GRID_SHARE: f64 = 0.10;

/// **32 tract sequences including the reference**, where the ordinary path allows six
/// (`candidate_alleles.md` §2.1). Settled by the owner 2026-08-24 on spec §12's Q2: a
/// repeat tract carries more real alleles than a SNP does, six was inherited from a
/// SNP/indel setting with nothing behind it for tracts, and HipSTR — the nearest
/// comparator — has no allele limit at all, admitting every sequence that clears a
/// per-sample test and abandoning the locus only if the haplotype product exceeds 1,000.
/// **Soft, and still unmeasured**: 32 is a judgement bounded by cost, since a locus over
/// `A` alleles has `A(A+1)/2` diploid genotypes and the loop scores every sample against
/// every one — 528 here against 21 at six, and 2,080 at 64.
pub const DEFAULT_MAX_CANDIDATE_ALLELES_SSR: MaxCandidateAlleles =
    MaxCandidateAlleles::new_or_panic(32);
```

### 2.3 What this path returns beyond the shared bundle

```rust
/// The shared selection, plus the two numbers the genotype prior needs and this path has
/// already computed.
pub struct SsrLocusSelection {
    pub selection: LocusSelection,
    /// **Parallel to `selection.alleles`**, reference included at index 0: each candidate's
    /// repeat count.
    ///
    /// `calling_priors.md` §5's `fill_ssr_seed` takes exactly this slice and **must not
    /// recompute it** — floor division and the motif have to land on one integer, and two
    /// producers is how they stop doing so (spec §3).
    pub repeat_counts: Vec<u32>,
    /// The ladder's mode, carried through for the same builder's `modal_repeat_count`.
    pub modal_repeat_count: u32,
}
```

---

## 3. The interface

```rust
/// Narrow one repeat tract's allele table to the tract sequences worth calling over
/// (spec §4, §5).
///
/// **Three passes over the locus, in this order**, because each needs the one before:
/// build the ladder and the per-sample length histograms; nominate rungs per sample and
/// union them; then admit the sequences on nominated rungs that clear the support rule, apply
/// the shared cap, and fill the leftover.
///
/// The reference tract is admitted first, exempt from the bar and the cap, and is admitted
/// **before** the periodicity test so a non-periodic locus still returns a usable table
/// (spec §7).
pub fn select_ssr(
    observation: &CohortObservation,
    detail: &SsrDetail,                    // §1 — no producer yet
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) -> SsrLocusSelection;
```

### 3.1 Nomination — the contract, per sample

For each covering sample, over its **spanning** reads only (partials are held apart and are read
by nobody here — [`candidate_alleles.md`](candidate_alleles.md) §3.1):

1. a repeat count is **nominated** when the sample's reads at it clear
   `config.shared.support` against that sample's spanning reads at the locus — the shared
   predicate, `MinAltReads::reached_by`, with the rung's read total as the numerator;
2. the top `ploidy` nominated counts by support are **promoted**;
3. if fewer than `ploidy` cleared it, each promoted count's `±1` neighbours are promoted too
   **where some sample's reads reached that length** — production's `occupied` test
   ([`candidate_set.rs:221`](../../../../src/ssr/cohort/candidate_set.rs)), kept so that nothing
   here invents a length;
4. the cohort's promoted set is the **union** across samples.

**Production's `is_clear_peak` is not called and is not ported**
([`rung_ladder.rs:274-288`](../../../../src/ssr/cohort/rung_ladder.rs)); spec §4.1 measures what
it costs. The `±1` rescue *is* ported, and it is the one part of production's nomination that
survives.

### 3.2 Admission within a promoted rung

Every sequence on a promoted rung faces the shared support rule, asked of the **sequence**: no
representative is privileged and no recurrence term applies (spec §5). That is the whole rule —
production's `cohort_alleles` and its three sibling constants have no counterpart here
([`candidate_set.rs:169-191`](../../../../src/ssr/cohort/candidate_set.rs)).

### 3.3 Periodicity

**Asked per sample, and one periodic sample is enough.** A sample is non-periodic when more than
`max_off_grid_share` of its spanning reads sit at a length that is not a whole number of motif
units from the ladder's mode; the locus is `SelectionVerdict::NotPeriodic` only when no sample is
periodic (spec §7). Production asks it of the cohort's reads pooled
([`candidate_set.rs:114-145`](../../../../src/ssr/cohort/candidate_set.rs)), and measures the
offset in **bases** while keying its ladder in **units** — ng measures both in units (spec §3).

A `NotPeriodic` locus returns the reference tract alone, that verdict, an empty leftover and a
`repeat_counts` of length one.

---

## 4. Design decisions — decided

- **The ladder is ported; the peak test is not.** — spec §4.1, measured: 33–78% against 97–100%
  on offering both alleles of an adjacent-length heterozygote, with fewer candidates.
- **The `±1` rescue and its `occupied` test are ported unchanged.** — spec §4.
- **One rule for rungs and for sequences within a rung.** The sibling bar of
  [`candidate_alleles.md`](candidate_alleles.md) §2.1 is asked twice with two numerators, not
  written twice. — spec §5, measured: 35% against 86% on the class that is 43% of HG002's
  heterozygous tracts.
- **The depth gate, the three sibling constants and the cap of 24 are dropped**, so
  `SsrSelectionConfig` has three fields and production's `CandidateCfg` has seven
  ([`candidate_set.rs:80-90`](../../../../src/ssr/cohort/candidate_set.rs)). — spec §6.
- **Periodicity is per sample.** — spec §7.
- **`repeat_counts` and `modal_repeat_count` are returned, not recomputed downstream.** One
  producer for an integer two modules must agree on. — this doc §2.3; the spec states the
  coupling, not the field.
- **`ploidy` comes from `FrozenParameters`**, never a constant in this module. Production hard-
  asserts diploid at its EM ([`em.rs:489`](../../../../src/ssr/cohort/em.rs)); ng does not, and a
  polyploid run changes only how many rungs a sample promotes.
- **The scratch is the shared `SelectionScratch`**, with the ladder's buffers added to it — one
  per-worker allocation for both paths. — [`candidate_alleles.md`](candidate_alleles.md) §2.4.

---

## 5. Reconciliation with existing code

| this doc's name | existing code | how it converges |
|---|---|---|
| `RepeatLadder` | [`Rungs` / `build_rungs`, `rung_ladder.rs:291`](../../../../src/ssr/cohort/rung_ladder.rs); the keying at [`:301-316`](../../../../src/ssr/cohort/rung_ladder.rs) | **shape ported**, keyed the same way. ng's holds table indices where production's holds cloned sequences ([`RungSeq`, `:40`](../../../../src/ssr/cohort/rung_ladder.rs)) — the merge already owns the bases |
| the per-sample length histogram | [`sample_histogram`, `rung_ladder.rs:262`](../../../../src/ssr/cohort/rung_ladder.rs) | same fold, over the merge's rows instead of `seq_counts` |
| nomination | [`assemble_candidates`, `candidate_set.rs:194`](../../../../src/ssr/cohort/candidate_set.rs) (top-`ploidy` at [`:231-238`](../../../../src/ssr/cohort/candidate_set.rs)) | top-`ploidy` and the union ported; the predicate replaced (§3.1) |
| the `±1` rescue and `occupied` | [`candidate_set.rs:239-258`](../../../../src/ssr/cohort/candidate_set.rs), [`:221`](../../../../src/ssr/cohort/candidate_set.rs) | **ported unchanged** |
| the clear-peak test | [`is_clear_peak`, `rung_ladder.rs:274`](../../../../src/ssr/cohort/rung_ladder.rs) | **not ported.** Kept only in the differential oracle's test arm (below) |
| same-length promotion | [`cohort_alleles`, `candidate_set.rs:169`](../../../../src/ssr/cohort/candidate_set.rs) | **not ported**; replaced by §3.2 |
| periodicity | [`is_periodic`, `candidate_set.rs:114`](../../../../src/ssr/cohort/candidate_set.rs) | shape ported, denominator moved per sample, offset measured in units |
| the reference seeded first | [`candidate_set.rs:200-202`](../../../../src/ssr/cohort/candidate_set.rs) | ported; it is [`candidate_alleles.md`](candidate_alleles.md) §1's invariant, held by `CandidateAlleles::new` |
| the depth gate, sibling constants, cap of 24 | [`CandidateCfg::dev_default`, `candidate_set.rs:80-90`](../../../../src/ssr/cohort/candidate_set.rs) | **not ported** (§4) |
| the motif and its period | [`Motif`, `types.rs:1107`](../../../../src/ng/types.rs), `period()` at [`:1138`](../../../../src/ng/types.rs) | reused; reached through `SsrDetail::motif` ([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)) |
| the prior's consumer | `fill_ssr_seed`, [`calling_priors.md`](calling_priors.md) §5 | `repeat_counts` and `modal_repeat_count` fill its first two arguments exactly |
| the junk term's denominator | `reachable_length_count`, [`read_likelihoods.md`](read_likelihoods.md) §4.1 | computed **from** the candidate set this returns; no field crosses, but widening the set changes it (spec §8) |

---

## 6. Open items

- `OPEN (blocked):` **`CohortObservation` carries no `LocusKind`** (§1). One field, the merge's.
- `OPEN:` **`DEFAULT_MAX_OFF_GRID_SHARE = 0.10` has never been measured by anyone** — spec Q1.
  Nothing in the code shape depends on the value; it is a config field with a named default.
- `OPEN:` **whether the shared cap binds at a tract in a large cohort** — spec Q2. A tract carries
  more real alleles than a SNP does, so this is where the shared document's extrapolation is most
  likely to come true. No type changes either way.
- *Impl-time:* whether the per-sample length histogram is a `BTreeMap` (production's) or a small
  dense vector over the locus's occupied range. The range is a handful of counts at almost every
  tract; measure before choosing.
- *Impl-time:* whether `RepeatLadder` is public or crate-private. It is written public here because
  the differential oracle's test arm needs to build one; if the test can go through `select_ssr`
  alone, make it private.

---

## Test & bench shape

Tests beside the file, per the repo rule; spec §13 lists what they assert. Two shapes belong to
this document rather than to the spec:

- **The differential oracle is a test-only arm, not a shipped config.** Spec §10 asks for
  production's three rules switchable so the differential has both ends. That switch lives in the
  test module — a `#[cfg(test)]` re-implementation of `is_clear_peak`, the cohort depth sum and
  the sibling bar, driving the same fold — so the shipping binary carries one rule and the
  comparison still has a failing state at both ends. Putting the old rules in
  `SsrSelectionConfig` would ship a configuration nobody should set.
- **The nomination test that matters most is the one production cannot pass:** a sample with 150
  reads at repeat count 10 and 150 at 11 nominates both. It asserts the difference from
  production, not a value, so it stays meaningful when the constants move.

The regression anchor is spec §13's last entry: HG002 at 30× through the whole caller with real
candidates, scoring at least the 98% adjacent-length and 86% same-length recall the spec measures.
The harness exists — production's `ssr-pileup` into `examples/ssr_slip_dump`, scored against
`benchmarks/ssr_hg002/truth/`.
