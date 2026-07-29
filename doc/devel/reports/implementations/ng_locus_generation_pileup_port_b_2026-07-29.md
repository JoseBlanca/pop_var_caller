# ng generic locus generator — the port, Milestone B: the differential, and proving it can fail

**Date:** 2026-07-29 · **Plan:**
[locus_generation_pileup_port.md](../../ng/impl_plan/locus_generation_pileup_port.md) steps B1–B3 ·
**Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md) §3, §8, §13.1

Implementation report for Milestone B of plan 2 of 3. Three commits, one per step — B2 unbundled
because "the fixture has been shown to discriminate" is a deliverable, not a side effect.

## 1. Plan

Prove that ng's copy of the walker **computes** what production's computes, where Milestone A proved
only that it **is** production's text. Then prove the proof can fail. Then run it at scale.

This is the last moment the baseline can be banked: plan 3's first commit makes the two walkers
differ on purpose, and the stage-1 differential dies with it.

## 2. Assumptions

One, and it was wrong in an interesting way: that a randomized differential over legal reads would
run to completion. It does not — production panics on ~4.2% of generated cases (§4). That changed the
harness's shape, not the plan's intent, and is recorded as a deviation in §6.

## 3. Changes made

### B1 — the stage-1 harness (`0a77b2e`)

[pileup/parity.rs](../../../../src/ng/locus_generation/pileup/parity.rs), `#[cfg(test)]`, in the
`delimit_parity` / `left_align_parity` shape.

- **One stream, fed to both.** Cases are generated as production's `PreparedRead` and converted to
  ng's through `from_production`, so both walkers see the same bytes in the same order *by
  construction*. Preparing separately would inject read preparation's uppercase divergence
  (`read_preparation.md` §6) into a comparison that is about the walk. Each walk builds its own
  `MockFasta` from the same bytes — the type is stateless, so that is equivalent to sharing one; the
  real-data test genuinely lends one `RefSeqFetcher`, whose accessor is not.
- **Byte-identity is well defined** because ng's copy still emits production's `PileupRecord`, whose
  hand-written `PartialEq` compares the two `f32`s by bits — so `finalise()`'s `NaN` placeholders
  compare equal and the comparison is total.
- **The generator** is randomized over four seeds (`SplitMix64`, written out rather than pulled from
  a crate, so a failure is reproducible from the source in front of you), with
  `PVC_PARITY_CASES` for a soak. It varies op mix, strand, quality, mapping quality, adaptor
  boundaries, mate pairing, contig, and the walker config — including tiny column caps, because the
  cap never fires at these depths otherwise.
- **`the_generator_exercises_what_the_port_can_break`** asserts the generator *reaches* each
  behaviour rather than assuming it. Measured from production's own `RunSummary` wherever it has a
  counter, so it reports what the walker did rather than what the generator intended.

### B2 — shown to fail (`f39f125`)

Five behaviours, mutated in ng's copy one at a time, each required to fail the differential and then
reverted. The table lives in `parity.rs`'s module doc so the exercise can be re-run from the source.

The table below is the **re-run** after the review changed the generator (§7); the original five all
failed too, at different case indices.

| # | behaviour | mutation | first divergence |
|---|---|---|---|
| 1 | mate-overlap reconciliation | early `return` from `resolve_mate_overlap_at_pos` | seed 0 case 2, item 6 |
| 2 | adaptor masking | `base_in_adaptor` always `false` | seed 0 case 0 — 27 records vs 25 |
| 3 | record widening | `widen` extends only `alleles[0]` | seed 0 case 18, item 17 |
| 4 | the subtract-then-add re-fold | the `subtract_contribution` half dropped | seed 0 case 0, item 16 |
| 5 | the column depth cap | `truncate(cap)` removed | seed 0 case 0, item 10 |
| 6 | **the panic *cause*** (added by the review) | ng's copy of the reachable `debug_assert!` replaced by an unrelated `panic!` | seed 0 case 15 |

All died inside the **first nineteen cases of the first seed**, against a default of 400 × 4.

**Mutation 5 is the one worth keeping.** It leaves `column_depth_truncations` incrementing, so a
differential that compared only the `RunSummary` — a plausible way to write this harness — would
have passed it. The records caught it.

**Mutation 6 exists because the review demonstrated it passing** against the first version of this
harness. See §7.

### B3 — at scale (`fcf0c3e`)

`ng_walks_identically_to_production_on_real_reads`: `#[ignore]`d, env-driven, ingesting real
alignments through ng's own step 1 and preparing them **once** with ng's `LeftAlignPreparer`.
Production's stream is that same stream through the new `PreparedRead::into_production`. One
`RefSeqFetcher` is lent to both walkers.

Re-run after the review, since the generator changed:

| run | scale | result |
|---|---|---|
| synthetic, release, `PVC_PARITY_CASES=5000` | 20,000 cases | **1,010,515 records, 0 divergences** |
| synthetic, debug, `PVC_PARITY_CASES=2500` | 10,000 cases | **492,224 records, 0 divergences**; 415 cases (4.2%) panicked in *both*, with the same message |
| GIAB HG002 10×, `chr1:1000000-1400000` | targeted TR bundle | **4,600 records, 0 divergences** |
| GIAB HG002 300×, `chr1:100000000-120000000` | 20 Mb | **137,591 records, 0 divergences** |
| tomato CRAM `SRR7279481.p1`, `SL4.0ch01:3406886-3506886` | 100 kb | **96,260 records, 0 divergences** |
| tomato CRAM `SRR7279481.p1`, `SL4.0ch01:13806669-15092603` | 1.3 Mb | **198,673 records, 0 divergences** |

437,124 records of real sequencing data across two organisms, a BAM and a CRAM, 10× and 300×.

## 4. The production defect this milestone found

**`apply_events_to_ref_into`'s `debug_assert!` is reachable on a legal read stream**, and the
release behaviour behind it is silently wrong.

Spec §8 already records that `events_overlapping` **does not clip a deletion to the window** — "one
anchored before the record can report an anchor below `record_pos`". What it does not record is that
production asserts the opposite. The mechanism, shrunk from a generated case to three reads and
pinned as `both_walkers_panic_on_a_deletion_anchored_before_its_record`:

1. a mate carries a deletion anchored at 17, spanning 18–22;
2. at 17 it overlaps its own mate, and mate-overlap reconciliation **in the indel regime collapses
   the pair to a single observation**, removing the contributor that carried the indel — so no
   record opens at 17;
3. another read's deletion opens a record at **19**, inside the footprint of a deletion that never
   opened one;
4. where the first mate matches again (23+) it folds into that record, and `events_overlapping`
   hands the fold its deletion anchored at **17**, two positions before `record_pos`.

Debug panics. **Release `saturating_sub`s the offset to 0**, applying the deletion's bases at the
record's first base — wrong allele bytes, no error, at exactly the long-deletion loci this port
exists to get right. ~4.2% of generated cases reach it.

Production is frozen, so this is **recorded, not fixed**. It is a Checkpoint B question (§8).

## 5. Validation

Container (`./scripts/dev.sh`), at each commit:

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics. *(It caught one real
  thing: `to_production` on a non-`Copy` type taking `self` by value. Renamed `into_production`,
  which also matches `AlignedRead::into_mapped_read`.)*
- `cargo test --all-targets --all-features` — **2648 passed / 0 failed / 5 ignored** after the review fixes (2644 before), plus the
  pre-existing `benches/psp_writer_perf.rs:386` panic.
- The soaks and real-data runs tabulated in §3.

`cargo doc --no-deps` not re-run: known-red on 11 unresolved intra-doc links. Both exceptions are
tracked under PROJECT_STATUS *Standing project-wide items*.

## 6. Deviations from the plan

1. **The differential represents panics rather than avoiding them.** The plan assumes both walkers
   run to completion; production panics on ~4.2% of legal generated inputs (§4). The alternative was
   to exclude that input class from the generator — which would have cut exactly the long-deletion,
   overlapping-mate inputs the port exists to get right, and hidden the finding. Instead
   `WalkOutcome` carries a panic flag and the differential requires the two walkers to agree on
   **which** inputs panic *and* on every record emitted before the panic. A verbatim copy must panic
   verbatim. This makes the parity claim stronger, not weaker.
2. **B3's real-data arm is an `#[ignore]`d in-crate test, not an `examples/` binary.** It needs
   `pub(crate)` conversions (`from_production`, `into_production`) that an example, being an external
   crate consumer, cannot reach. Env-driven so one test serves both organisms.
3. **B3's real-data arm wants `--release`.** Not an optimisation: in debug it would hit §4's
   assertion constantly on ordinary paired-end data and measure almost nothing. Documented on the
   test.
4. **`PreparedRead::into_production` is new**, `#[cfg(test)]`, and exists for exactly one caller —
   B3, which must hand one prepared stream to both walkers. Destructured, so a field added to ng's
   type stops it compiling rather than quietly ceasing to be carried.

## 7. Review, and what it changed

Reviewed the same day over the milestone diff:
[ng_locus_generation_pileup_port_b_2026-07-29.md](../reviews/ng_locus_generation_pileup_port_b_2026-07-29.md)
— 6 categories, **0 Blockers, 5 Majors, 10 Minors**, Approve-with-changes. All applied.

**Two reviewers went past reading, and that is where the value was.** The reliability pass re-ran two
of B2's five mutations and reproduced the table to the seed, the case index and the stream item — then
ran a **sixth mutation the table did not cover, and the whole suite stayed green.** The generator pass
re-implemented the generator in a scratchpad, verified fidelity against the harness's own printed
counter, and measured the population over 18,077 reads.

**The five Majors:**

- **The panic *cause* was discarded.** `WalkOutcome.panicked` was a `bool`, so replacing ng's copy of
  production's reachable `debug_assert!` with a semantically unrelated `panic!` passed —
  including the test whose own doc claimed to check ng "reach[es] the same precondition".
  `open_record.rs` carries eight distinct `debug_assert!`s. **This is the same pattern as the other
  four milestones, and this time it was in a claim I wrote into the module doc.** Now the *message*
  is captured and compared, the comparison runs **first** so a panic divergence is diagnosed as one,
  and the pinned defect is pinned by cause.
- **The generator emitted 5 of 9 `CigarOp` variants** — `=`/`X`/`H`/`P` never appeared, 0 of 18,077
  reads. `=`/`X` are minimap2 `--eqx` and DRAGEN output and share the `Match` arm at four cursor
  sites.
- **The real-data test had no record floor**: a walk dying on read one prints "1 records compared,
  zero divergences" and passes.
- **`reads_with_live_adaptor_boundary` was exactly `is_some()`** — measured 4542 = 4542, of which 126
  (2.8%) silence nothing. The one assertion written to prevent "a test that cannot fail" was one.
- **The `Err` half of the stream had no inputs behind it** — 0 walker errors over 1600 cases, while
  spec §3 claims the two `Result` streams are compared element for element.

Ten Minors applied, the substantive ones being: a named `SummaryCounters` built by two exhaustive
destructures (the old macro read fields *by name*, so a ninth `RunSummary` field would have silently
left the parity claim); reads placed at the contig end, where four bounds guards had been provably
inert and no record ever widened to a contig's final base; the mate-eviction window drawn small
enough to bite, and asserted (61 → 221); one generic `drive` replacing four copies of the
catch-unwind block, two of which had already drifted; and `into_production` given the runtime
coverage it had none of.

## 8. Checkpoint B — questions for the owner

1. **The production defect (§4).** Recorded, not fixed, because production is frozen. Worth knowing:
   plan 3 replaces the haplotype builder that contains it, so the *ng* side stops being wrong as a
   side effect — but production keeps mis-placing deletion bases at these loci until someone
   decides to touch it. Does it want a research note of its own, or an entry in the production
   defect list beside the partial-coverage ref-fill finding?
2. **Milestone A's four questions are still open** — the copy banners, the "46 tests" in the plan and
   the spec, the arch doc's file inventory, and whether `RefSeqFetcher` should be renamed and moved.
   None blocks plan 3.

## 9. What Checkpoint B hands to plan 3

- **The baseline is banked.** 1.4 M synthetic records and 437 k real ones, zero divergences, over a
  differential that has been shown to fail five ways. Plan 3's first commit is where the two walkers
  start to differ on purpose.
- **`parity.rs` dies with that commit, and it should.** What survives is narrower and is plan 3's to
  build: on loci where every folded read witnessed the whole footprint, the two must agree forever.
- **§4 is a gift to plan 3's `read_coverage` work**, not just a defect report: it is a concrete,
  reproducible case where production's fold reaches for bases a read did not witness, with the
  mechanism traced.
