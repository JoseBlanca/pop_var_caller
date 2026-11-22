# Refinement's memory, second cut: the same output in less (2026-11-22)

Branch `refine-memory-2`, from main c60a9bf6.

## Why

The second tomato25 run (2,176 samples) stopped in the writing pass at chromosome 1, 42.25 Mb: the
smallest window, read only part-way, already counted 93.3 GB against a 90 GB budget. Cut from the
run's kept spill and measured record by record, the whole window at the floor (0.5-Mb core and
margins) holds 246,848 records — 160 a kilobase, 43% of them with more than two alleles — and
counts about 123 GB; real memory ran about 25% above the count on top of that. Every cut here
changes no number the pass computes.

## Changes made

- [`records.rs`](../../../../src/neighbour_evidence/records.rs): **a record's likelihoods are its
  `PL` above the sample's lowest, in a `u16`**, read as `10^(−PL/10)` through a table of every
  `u16` made by the expression the `f32` was made by (`likelihood_of`): two bytes a genotype for
  four, the same `f32` read back. Past about 455 the `f32` is zero, so clamping at `u16::MAX`
  changes nothing.
- [`phase.rs`](../../../../src/neighbour_evidence/phase.rs): **each allele pair's path bits packed**,
  the paths' `towards_lower` bits then their `leans` bits, `⌈2·paths/8⌉` bytes — 9 for the 36 paths
  of a run's window, where two `u64`s took 16. The window's path count sets the packing
  (`Phasing::paths`); the budget counts the run's (`refine::PHASED_PATHS`).
- [`fragments.rs`](../../../../src/read_phase/fragments.rs): **read alleles in three bytes** where
  every allele of the record is at most 3 and every Phred at most 63 (`TightReadAllele`: the
  `u16` offset, then allele and Phred in one byte); four or seven bytes otherwise, as before.
- [`stage.rs`](../../../../src/neighbour_evidence/stage.rs), [`writer.rs`](../../../../src/vcf/writer.rs):
  **a core's rewritten lines go to the file 1,024 at a time** (`push_writing`, `finish_writing`,
  `StageStop`), where every one was held until the stage returned — about 9 GB in the failing
  window. The budget no longer counts a core's lines. One `DEBUG` line a window in the run's log:
  where, how many records, the bytes counted, the seconds taken.
- [`model.rs`](../../../../src/neighbour_evidence/model.rs): **a walk's blurred shares and likelihood
  matrix are computed where they are read** (`SiteTerms::fill`) into a buffer it reuses, rather
  than kept for every site of the window on every thread — a few hundred bytes to over a kilobyte
  a site a thread; and **the drawn definite alleles are a byte each** (`WalkSite::emission_shares`),
  where they were `f64` zeros and ones.
- [`refine_replay.rs`](../../../../examples/refine_replay.rs),
  [`spill_slice.rs`](../../../../examples/spill_slice.rs) (new examples): cut a region out of a
  kept spill, and replay refinement over it at the run's scale — its output, each window's peak
  memory and time, or each record's counted bytes; `REPLAY_FIXED` fixes the windows so builds that
  count differently solve the same ones.

## Tests added or changed

- `records::tests::a_likelihood_held_in_two_bytes_is_the_f32_it_stands_for`: every `u16`, bit for bit.
- `fragments::tests::read_alleles_are_held_in_three_bytes_where_their_alleles_and_phreds_fit`.
- The phase tests make their phase for 64 paths; the emission-shares test gives bytes.
- The recorded digests of the posteriors and the phase did not move.

## Validation

- `cargo test --release --lib`: all pass but the two tests that expect a debug-only panic, as on
  main; `cargo test --release --test '*'`: 33 passed.
- **The 2,176 tomato25 samples, chromosome 1 at 44.0–45.5 Mb from the run's spill**, every window
  0.5-Mb core and margins, 48 threads, main against this branch:

| | main | this branch |
|---|---|---|
| Output (11.6 GB, 130,383 records) | | byte-identical |
| Peak memory, first window | 50.4 GB | 39.4 GB |
| Peak memory, the run | 62.5 GB | 51.6 GB |
| Wall time | 1,705 s | 1,174 s |
| CPU time | 17.4 h | 10.7 h |

The walks were faster for holding less: their terms no longer travel through memory.

**Counted bytes** of the failing floor window (246,848 records): about 123 GB on main, 98.5 GB
here — a sample a record, 76.5 bytes of solve state held over the rounds (37 at a biallelic record,
259 at a six-allele one), 46 of read alleles, 28.4 of the core's sums and results, 24.3 of the
record's likelihoods and shares, 8.1 of compressed lines. The densest floor window between 40.5 and
46 Mb counts 101 GB. Real memory ran 4.6 to 10.6 GB above the count here (39.4 GB against 34.8,
51.6 against 41.0).

## Tradeoffs and follow-ups

- **The densest windows still do not fit 125 GB at the floor.** What is left is used as held; the
  next lever is fewer records a window there (owner's call, 2026-11-22).
- A dense index of each record's samples with read alleles would save about 2 bytes a sample a
  record: not done.
