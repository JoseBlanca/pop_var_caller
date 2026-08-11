# ng step 4, generic path — Milestone C: one locus → one cell

**Date:** 2026-08-06. **Branch:** `ng-parameter-estimation`.
**Plan:** `doc/devel/ng/impl_plan/parameter_prepass_generic.md`, Milestone C (C1–C3),
stopping at Checkpoint C. **Design:** `arch/parameter_prepass_generic.md` §2.2, §2.3, §3;
`spec/parameter_prepass_generic.md` §1, §4, §9, §12; research note §4.6.

## 1. What was built

The path from one locus to one cell, and the accumulator a sample's loci are poured into.

| commit | step | what landed |
|---|---|---|
| `1977c465` | C1 | `count_whole_site`, `count_by_read_group` — the only place that decides what an alternative read is |
| `4b01681d` | C2 | the depth cap: `CountedSite`, `hypergeometric_draw`, `seed_at` — **own commit, do not bundle** |
| `eccccf7c` | — | C2 amendment: `count_whole_site_by_library`, `SelectionWalk` (forced by C3) |
| `9d3e090d` | C3 | `GenericAccumulators`, `add_locus`, `AccumulationCounts`, `PloidyMap`, `InbreedingMode`, `WindowKey` |
| `1c7b4490` | — | the three-agent review applied |

New file: `src/ng/parameter_estimation/generic/accumulators.rs`. `depth_and_alt_reads.rs`
went from a doc stub to C1, C2 and the amendment.

## 2. The three decisions that carry the milestone

**What counts as an alternative read** — byte equality against the locus's own
`reference_bases`, over complete witnesses only. It is the per-read form of production's
`allele_index == 0` and needs no special cases: `bases` is in *read* coordinates, so an
insertion is longer than the reference slice and a deletion shorter, and neither compares
equal, which is right for both.

**One dependency is invisible and would have been catastrophic.** ng has two reference
readers: `RefSeq::fetch_into` returns canonical `{A,C,G,T,N}`, and `RawChromReader`
preserves soft-masked lowercase verbatim for the typed-region catalog's byte oracle. The
generic locus generator fetches through the first — checked before the comparison was
written. Had it used the second, byte equality would call **every read at a soft-masked
position alternative**: about half the human genome, and the repeat-rich half.
`byte_equality_is_the_rule_and_a_lowercase_reference_would_break_it` shows both sides.

**Subsample, do not rescale.** A site deeper than the ladder's cap keeps 124 of its reads
and counts the alternative ones among them, seeded from the locus position. Both
shortcuts were written and run against the oracle, and **the variance is what catches
them, by two orders of magnitude**: round-to-nearest gives every site the same answer
(variance 0), the stochastic round gives `f·(1−f)` ≈ 0.19, the truth is 23.4. All three
have the mean right to well inside the test's tolerance, so a mean-only assertion would
have accepted either.

## 3. What C3 forced on C2, and why it could not have been foreseen

`add_attributed_site` asserts a site's library attribution sums to that entry's own
alternative count. Above the cap, `count_by_read_group`'s per-group counts cannot serve:
each is its own subsample of its own group's reads drawn against its own depth, so their
alternative counts sum to something else. The assertion would have fired on the first
multi-library sample above ~124×.

Nothing in Milestone B or in C1/C2 could have caught it — each function was right about
its own question, and only the caller that uses both at once reveals the answers do not
compose. `count_whole_site_by_library` answers the windowed table's question from **one**
resumable `SelectionWalk`, so the per-group counts sum to exactly what a single call would
return, by construction rather than by two implementations agreeing.

**The two tables draw differently above the cap, deliberately.** They ask different
questions of the same site — one about each library's chemistry, one about the
individual's genotype. What matters is that each is internally consistent and that at one
read group they coincide exactly, which is what spec §12.6 needs.

## 4. Deviations recorded

- `count_*` take `&DepthBinEdges`, not a bare cap, so the cap and the ladder cannot
  disagree. No `max_site_depth` alias: `edges.max_depth()` *is* that function, and a
  second name is a second place to disagree.
- `CountedSite` replaces a bare `DepthAndAltReads` return, carrying the pre-cap depth so
  C3 can tally it.
- `count_whole_site_by_library` has no counterpart in the architecture (§3 above).
- `adjustments()` returns owned, where arch §3 declares `-> &AccumulationCounts`: the
  seam count is *derived* at read time by sorting, which is what keeps merge
  order-independent.
- `WindowKey` is a two-arm enum where the arch sketches a plain tuple key, so
  `InbreedingMode::Supplied` collapses the object honestly rather than through a sentinel
  key that lies.
- `AccumulationCounts` gains a fifth field, `shard_spans_overlapping`. The architecture
  implies one counter that "still merges"; counting the *loci* that overlap across a seam
  needs the loci, and a shard keeps only its span. Two counters say what is known.
- **The 1-based window arithmetic is settled here**: `(start − 1) / 100_000`, so every
  window holds exactly 100,000 positions. The naive division would leave the first window
  of every contig one base short, at a different resolution from all the rest. A locus is
  filed by where it *starts*.

## 5. What the reviews found

Three agents, eleven categories, ~60 mutations, **20 survived**. Three Blockers, and all
three were tests that could not fail. Full detail in the review synthesis; the two worth
carrying:

- **`merge` could drop the read-group table and all four counter sums, green.** A
  region-sharded walk would produce an empty read-group table — every error rate this step
  fits — with nothing to show it. Found independently by two agents.
- **`SelectionWalk`'s all-kept fast path was never resumed.** Dropping its population
  decrement leaves later groups drawing against a population one too high, and nothing
  panics, because the entry's attribution still sums to its own total.

Two findings were upstream or documentary and are corrected: arch §2.2's worked example
for round-to-nearest was on the wrong side of its own sign reversal (`round(1 × 124/500)
= 0`; the reversal is at depth 248), and "two groups is every sample with more than one
library" is false — 20 tomato samples carry three and four carry 7, 16, 16 and 42, which
was costing two allocations per locus from three groups upward.

## 6. Deferred, recorded rather than acted on

- **`AccumulationCounts` is one type with two meanings.** The stored value's
  `shard_spans_overlapping` is always zero; only `adjustments()`'s is real. That is what
  forces `merge` to enumerate the other four fields by hand (now exhaustively
  destructured, so a fifth cannot be forgotten). A `StoredCounts` / `ReportedCounts` split
  would say it in the type.
- **`CountedSite` means three slightly different things** across its three producers: on
  the by-group path it is a *group's slice* of a site, so summing `subsampled_from()`
  there would count groups under a name that says sites. No consumer does today.
- **`generic/mod.rs` now does four jobs** and splits `WindowIndex` from `WindowKey` across
  files; it still carries a note addressed to "whoever writes the division (Milestone
  C)", which is now written.
- **`add_locus` after `merge` is silently unsound** — `previous_end` is per-shard by
  design, so a locus added after a merge is checked against the wrong reach. No caller
  does this; the fold shape makes it unnatural.
- **The split-borrow alternative to `mem::take`** was built by a reviewer and passes: it
  removes both takes and the per-locus `Arc::clone` for 1–2 ns a locus. Worth taking when
  `add_locus` is next touched.

## 7. Validation

All via `./scripts/dev.sh`, at every commit:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --all-features` — **2,996 → 3,030 passed**, 0 failed,
  5 ignored. `ng::parameter_estimation` holds **128** tests, 94 → 128 over the milestone.
- `cargo doc --no-deps --lib` — **12 unresolved links, the pre-existing baseline**. It
  caught one regression the other three gates did not: the C2 amendment linked
  `DepthAltHistogram::add_attributed_site` from a module that does not import it.

`cargo test --all-targets` remains red through `benches/psp_writer_perf.rs:386`, in frozen
`src/psp/` code byte-identical to the branch point. The owner has ruled it out of scope.
