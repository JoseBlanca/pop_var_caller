# Fixes applied — the cohort merge's "which sample comes next"

*2026-08-18. Against [the review](ng_cohort_merge_walk_order_2026-08-18.md).*

**Everything actionable was applied, including the one that replaced the structure.** Three
Majors, five Minors, the nits.

- **M1 — the heap became a tournament tree over the covering samples**, which is what the
  review measured as fastest at every cohort size from 10 up and at both the dense and the
  sparse end. Re-measured here after the port, median of seven repeats: **15.9 µs at 63 samples
  against the heap's 19.2, 81.9 µs at 250 against 119, 391 µs at 1,000 against 773, and 2.31 ms
  at 3,000 against 3.10** — and against the scan this work started from, 3.6× to 44×.
  **One hardening on top of the reviewer's patch:** exhaustion is now *counted* rather than read
  off the sentinel key that sorts a spent leaf last, so a real observation on contig `u32::MAX`
  at position `u64::MAX` cannot be mistaken for an exhausted sample.
- **M2 — the hazard is deleted rather than asserted.** `take_head` takes no sample: it reads the
  winner out of the structure, so consuming a sample other than the one shown is
  unrepresentable rather than `debug_assert!`ed.
- **M3 — the oracle's generator now sweeps the committed range**: 1 to 400 samples, a sample
  allowed to cover nothing, and first bases spread over the whole stretch so one sample's
  observations can lie wholly before another's. Three contigs rather than two.
- **Mi1** — `the_tournament_names_the_lowest_sample_index_when_two_samples_start_together`
  pins the tie-break directly on the private helper. It dies under a flipped key, which no
  other test in the file does.
- **Mi2** — `the_tournament_holds_one_live_leaf_per_unspent_sample` asserts the structure the
  cost rests on: one live leaf per unspent sample at every step of a 64-sample walk, and every
  observation consumed exactly once. A scan-based rewrite closes the same loci and cannot
  satisfy this.
- **Mi3** — the key is a named 16-byte `HeadKey` laid out `position, contig, sample`, with `Ord`
  restoring genome order; the doc says the field order is layout and the `Ord` impl is meaning.
- **Mi4** — the probe reports the median and the min–max of seven repeats, raises the walk count
  so each timed block is long beside the jitter, and prints one column instead of two, saying
  why the second is gone.
- **Mi5** — the wrong ratio is corrected: the doc now gives all four cohort sizes and says the
  4×/14× is the step between adjacent ones.
- **The identity check** the reliability review asked for is folded into the randomised oracle:
  every sample is handed back its own observations, once each, in its own order, none outside
  its locus's ground.

## Not done, and why

- **A criterion bench under `benches/`.** The repo has ten and the merge walk is not among them,
  so nothing in CI notices this cost regressing. What guards it today is the structural test
  above. Recorded as open rather than added, because a bench is a project-level decision about
  what CI runs.

## Mutation re-run

Ten mutations against the tournament, **all ten killed**: the match keeping the winner instead
of the loser; a key ordering position before contig; the tie-break flipped; a spent leaf sorting
first; the seed above the real keys; no rebuild pass at construction; the replay starting a
level too high; a spent leaf not counted; the wrong sample's leaf; and empty samples given a
leaf. Nine earlier mutations against the heap were also all killed.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` → `168 passed; 0 failed`; the whole library suite green.
