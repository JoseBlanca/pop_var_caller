# ng cohort merge — which sample comes next: a tournament, not a scan over the cohort

*Implementation report, 2026-08-18. Not a step of
[the plan](../../ng/impl_plan/cohort_merge.md) — the owner asked for it at Checkpoint D,
after the measurement below was raised.*

## 1. What was wrong

`LocusCloser` walks k samples' observations merged by position, and to do that it must
answer, for every observation, **"which sample comes next?"**. It answered by scanning all k
samples' current positions, so the walk grew with the *square* of the cohort. Measured in a
release build on one 20-base building region with every sample carrying a record at every
position of it — which is what the generic mint produces wherever a sample has reads:

| cohort | scanning the cohort |
|---|---|
| 63 samples | 57 µs |
| 250 | 793 µs |
| 1,000 | 11.6 ms |
| 3,000 | 101 ms |

Four times the cohort for about fourteen times the time. At 63 samples that is invisible,
which is why it had never shown; at the thousand-to-three-thousand end the caller commits to,
it decided whether the merge ran at all — about **57 days of walking for an 800 Mb genome on
one thread**, which sixteen threads would only bring to three and a half days.

## 2. What it is now

A **tournament tree** — one leaf per sample that has anything in the walk, every internal node
holding the loser of the match below it, the winner at the root. Taking an observation writes
that sample's next beginning into its leaf and replays the matches up to the root, **one
comparison a level**.

**A binary heap was built and measured first**, and the tournament beat it, so the heap is
recorded here rather than shipped. Three things it was paying for that the merge does not
need: a 24-byte key where 16 will do (`GenomePosition` puts a `u32` before a `u64`, so a third
of it is padding), a pop followed by a push where one replace-top would do, and a sift-down
that compares *both* children at each level to decide where to sink. The tournament compares
once a level, because the tree already remembers who lost — and it can, because **k is fixed
for the whole walk**: which samples have observations is settled at construction and no sample
ever joins.

**The leaves are the covering samples, not the whole cohort**, and that is what makes it win
at the sparse end too. Over the whole cohort it is *slower* than a heap where most of the
cohort has nothing in the region — 74 µs against 38 at 3,000 samples with 1% covering, because
it builds 3,000 leaves to merge 30 — and keyed to the covering samples it is the fastest of
everything the review measured there, 28 µs.

## 3. What it costs now

The same regions, median of seven repeats, release build
(`examples/ng_cohort_merge_walk_cost.rs`):

| cohort | scan | heap | tournament | against the scan |
|---|---|---|---|---|
| 1 sample | 0.40 µs | 0.45 µs | 0.35 µs | — |
| 10 | 2.57 µs | 2.37 µs | 1.92 µs | 1.3× |
| 63 | 56.9 µs | 19.2 µs | 15.9 µs | **3.6×** |
| 250 | 793 µs | 119 µs | 81.9 µs | **9.7×** |
| 1,000 | 11.6 ms | 773 µs | 391 µs | **30×** |
| 3,000 | 101 ms | 3.10 ms | 2.31 ms | **44×** |

The cost now grows about as fast as the cohort itself rather than as its square. For an 800 Mb
genome at 3,000 samples that is roughly **25 hours of walking on one thread against 57 days**.

**Nothing else in the merge changed**, and the ordering is most of what is left to change: the
review measured an oracle that already knows the merge order and pays nothing to find the next
sample, and it runs the same fixture in 7.3 µs at 63 samples and 1.03 ms at 3,000 — so
ordering was about two thirds of the walk before this, and the tournament closes about half of
that gap.

## 4. The output does not change, and that is what most of the work went on

- **The 168 tests of the module pass unchanged**, and this change was made to pass them
  without touching one.
- **A new randomised differential** compares the walk against an oracle written a different
  way — every observation in one list, sorted, then chained, with no merge in it at all — over
  300 random cohorts. Widened after the review found its first generator reached none of the
  shapes that separate a tournament from a scan: it now sweeps 1 to 400 samples, lets a sample
  cover nothing, and lets one sample's observations lie wholly before another's.
- **It also checks by identity** that every sample was handed back its own observations, once
  each, in its own order, none of them outside its locus's ground. Per-sample *counts* cannot
  see a member window shifted by the same amount at both ends.
- **Nineteen mutations** were run against the two versions and all nineteen fail at least one
  test: ten against the tournament (the match keeping the winner instead of the loser, a key
  ordering position before contig, the tie-break flipped, a spent leaf sorting first, the seed
  above the real keys, no rebuild pass, the replay starting a level too high, a spent leaf not
  counted, the wrong sample's leaf, and giving empty samples a leaf), and nine against the heap
  before it.

## 5. Three things the review changed beyond the structure

- **The hazard was deleted rather than asserted.** `take_head` used to take the sample as an
  argument and `debug_assert!` that it was the one the caller had been shown. The reviewer made
  them disagree and measured what follows in a *release* build: one wrong locus is emitted —
  holding one sample's record while another's two observations vanish from the run — before the
  walk dies on a key with no head. It now reads the sample out of the structure, so the
  disagreement is unrepresentable.
- **Exhaustion is counted, not inferred.** A spent leaf holds a sentinel key that sorts last.
  Reading "the walk is finished" off that sentinel would mistake a real observation on contig
  `u32::MAX` at position `u64::MAX` for an exhausted sample; a counter of spent leaves costs
  one increment per sample per walk and removes the question.
- **The probe reports a median and a spread** over seven repeats, and one column instead of
  two. Repeated runs of one unchanged binary at 3,000 samples gave 3,342, 3,412 and 4,438 µs —
  a 30% swing a single mean cannot tell from a code change; and the second column (a region
  with one substitution in it) measured within that swing of the first at every cohort size, so
  it implied a difference the measurement does not carry.

## 6. What is still open

- **No `benches/` guard.** The repo runs ten criterion benches and the merge walk is not among
  them, so nothing in CI would notice this cost going back up. What guards it today is
  structural: a test asserts one live leaf per unspent sample and one consumption per
  observation, which a scan-based rewrite cannot satisfy.
- **Construction is 1% of the walk** at 63 samples and 0.35% at 3,000 (measured), so making the
  closer reusable across building regions would buy at most that — except on a region almost
  nobody covers, where it is 7%.

## 7. Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` → `168 passed; 0 failed`; the whole library suite
green.
