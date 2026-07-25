# ng — read preparation (STR path) — RETIRED: there is no STR read-preparation step

*Retired 2026-07-25. **The STR path has no read preparation.** Read preparation is a per-read,
*locus-independent* transform that canonicalises a read against the reference around its own span
(left-align / BAQ / re-align → `PreparedRead`) and is reused across every locus the read overlaps — see
[`read_preparation.md`](read_preparation.md) §1. The STR per-read operation is not that: it aligns a
read **against a specific tract** to read out what the read shows **about that locus**, so it is
observation generation, not preparation. What this file used to spec now lives in its true homes, and
this redirect keeps older section links resolvable.*

## Where this file's content lives now

| former section here | what it covered | live home |
|---|---|---|
| §1, §3 | the observation shape; **partial (lower-bound) observations** — reads that anchor one flank and run off mid-tract, kept as `Partial` evidence (production's `BorderOffEnd`, no longer dropped) | [`locus_generation_ssr.md`](locus_generation_ssr.md) — the `SampleLocusObservations` / `ReadCoverage::{Complete, PartialLeft, PartialRight}` design, **built** (`src/ng/locus_generation/`) |
| §2 | the per-read tract alignment (the delimiter) — the "ruler, not a scorer" Viterbi/best-path alignment | [`alignment.md`](alignment.md) §4.2 — the repeat-aware best-path aligner, **built** (`src/ng/alignment/ssr_best_path_flat_gap.rs`, algorithm 4). The STR generator **calls** it, per locus, per read |
| §4 | the no-observation reasons (no border anchored, low quality, window truncated) — tallied, reported | [`locus_generation_ssr.md`](locus_generation_ssr.md) |
| §6 | the **censored likelihood** a `Partial` needs before step 7 can consume it (hard inequality, `1/allele` dilution, conservative discount — GangSTR's `FlankingClass`) | step 7 (`ReadLikelihoodModel`); flagged in [`locus_generation_ssr.md`](locus_generation_ssr.md) |
| §8 | open questions — where alignment is invoked (per locus, per read, by the generator); widening; whether partials pay | closed / owned by [`locus_generation_ssr.md`](locus_generation_ssr.md) |

**One consequence recorded here because it is easy to miss:** the read-selection gate upstream must
admit partially-covering reads, or the `Partial` observation class is unreachable. That is
[`locus_generation_ssr.md`](locus_generation_ssr.md)'s concern.
