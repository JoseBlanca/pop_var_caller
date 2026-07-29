# Code review — ng generic locus generator, Milestone C (C1–C4)

**Date:** 2026-07-29 · **Scope:** `fb8dde0..52d99c5` on `ng-generic` (C1 `a5c0203`,
C2 `d23e2b0`, C3 `c89d596`, C4 `52d99c5`) · **Fixes:** `94758d7` ·
**Impl report:** [ng_locus_generation_pileup_generator_c_2026-07-29.md](../implementations/ng_locus_generation_pileup_generator_c_2026-07-29.md)

Five category agents, **each in its own git worktree** — the isolation the
fan-out skills started mandating this session. Zero collisions again; every
result first-hand, nothing needing serial re-verification.

| agent | categories |
|---|---|
| reliability | `reliability` + the challenge-tests pass |
| errors + defaults | `errors`, `defaults` |
| idiomatic + smells + sharing | `idiomatic`, `smells`, `unsafe_concurrency` |
| structure + naming | `refactor_safety`, `module_structure`, `naming` |
| intent + cost | `extras` — spec conformance, the memory property, the cost of a region |

**Verdict: 2 Blockers, 6 Majors, and a long tail of Minors.** Both Blockers were
found independently by more than one agent. Everything applied is in `94758d7`;
what was carried is in the impl report's Checkpoint C section.

## The two Blockers

**1. A shed error outlived the region that shed it.** Found by the errors agent
and, separately, by the structure agent. The stream hands the walker an
infallible item type, so its fatal errors are shed into a cell and collected a
call later — but the cell lived on the generator-lifetime `ReadPreparation` and
was read on exactly one of three exit paths. Reproduced twice: a chr2 read's
preparation failure reported against a healthy **chr1** region after it had
emitted all thirty of its loci; and, with no region following, an error still
sitting in the cell when the generator dropped. Both are silent failures in the
strict sense — one misattributes, the other discards.

**2. Neither `Ok(None)` nor `Err` was terminal.** Found by three agents and by
the orchestrator. `end_walk` cleared the walk and not the region, so the call
after either re-opened the query and re-walked: `reads_admitted` 1 → 2 and
`chain_allocations` 1 → 2 — **the same fragment with two chain ids**, which is
precisely what a run-lifetime allocator exists to prevent. Latent through
`GeneratorSet` (which `take()`s its current region) and live through the
generator's own `pub fn next_locus`, which every one of the module's tests
drives. The reliability agent noted that the sibling `NoLoci` *is* idempotent
past the end and the public wrapper declares `FusedIterator` — so two impls of
one trait disagreed.

## The Majors

1. **Four knobs had an unmeasured envelope.** Zero is accepted by `check()` and
   then means three different failures, the worst being `max_snp_column_depth:
   0` — walks "successfully", emits **zero loci for a covered region**, and
   counts the 32 truncations that explain it into a struct the dispatcher cannot
   read. *Fixed:* every knob has a floor.
2. **The halo width was unpinned** — a halo of `max_record_span / 2` passed all
   25 tests. *Fixed:* the far read now sits at the halo's far end.
3. **The stop rule's comparison was unpinned** — `>` → `>=` passed 188 tests and
   drops the locus anchored exactly at `region.end`: a one-base hole at every
   region boundary. *Fixed*, and the first replacement fixture **also failed to
   catch it** (the rule's open-record half covers for the mutation whenever
   coverage runs on into the boundary); the shipped one starts coverage *at* the
   boundary.
4. **Two of `fold_region_walk`'s three rules were unpinned**, on the two counters
   no BAM fixture can reach. *Fixed:* the fold is exercised directly.
5. **`LocusGenerationError::Reads` has two origins and denied one of them** —
   C2's `open_walk` made the open failable, and the variant's doc said the open
   had already succeeded. *Doc fixed; the split carried.*
6. **The read-query accessor factory is called once per file per region** — the
   shape of spec §8's ~564k-opens trap, in the one accessor the generator cannot
   hold for the run. No non-test caller exists yet. *Documented on the field;
   carried.*

## What the intent-and-cost agent measured

- **Conformance holds at scale.** 199,990 loci, **identical** across 1, 2, 10,
  50, 200 and 500 regions, at both halo widths and with the stop rule ablated —
  the tiling property, tested rather than argued.
- **The memory property holds, and nothing was watching it.** Instrumented pull
  counts show the query consumed one read at a time; then an ablation collecting
  it into a `Vec` left **the whole library green, parity included**. That gap is
  now closed by a test.
- **A region costs ~0.12 ms, flat from k=10 to k=500**: `T(k) ≈ 81.5 ms + 0.12 ms
  × k` over a 200 kb contig at depth 3. In transferable terms, **a region costs
  about 290 loci of walking to set up**, so regions under ~300 bp cost more to
  open than to walk. Peak bytes are **flat** at 257,621 whether a region is 400 bp
  or 200 kb — the depth-shaped footprint, measured.
- **The stop rule is worth 19× the per-region constant and 8.7× total** (1,242 ms
  → 142 ms at 500 regions) **while changing not one emitted locus.** The halo is
  free because the stop is there.

## Declined

**Renaming `generator.rs` to `region_walk.rs`**, proposed on the rule that
renamed `driver.rs` → `genome_walk.rs`. The file is named for the type it
defines, as `ssr.rs` is; a second `*_walk.rs` beside `genome_walk.rs` would be
the less clear of the two.

## The pattern that keeps recurring

Milestone C shipped **two** tests that could not fail — C2's original stop test,
and the replacement written for the review's own finding. Both were caught only
by mutating the code the test named. The count for the branch is now nine.
