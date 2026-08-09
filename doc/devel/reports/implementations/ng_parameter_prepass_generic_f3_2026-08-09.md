# ng step 4, the SNP/indel path — F3: the identities on real alignments

**Date:** 2026-08-09. **Plan:** F3. **Design:** `arch/parameter_prepass_generic.md` §9,
`spec/parameter_prepass_generic.md` §12.6.

Tests only — no production behaviour changed. This is the first time step 4 reads a real
alignment file.

## What was built

`generic/real_alignments.rs`, a `#[cfg(test)]` module of four `#[ignore]`d tests driven by
`PVC_PREPASS_FASTA`, `PVC_PREPASS_READS` and `PVC_PREPASS_BED`, following
`locus_generation/pileup/parity.rs`'s convention (environment, `--release`, invocation in the
doc comment). Each test builds the real stream — typed-region catalog over the BED,
`PileupGenerator` in a `GeneratorSet`, `SampleReads`, `LeftAlignPreparer` — and feeds
`GenericAccumulators::add_locus`.

| test | what it asserts |
|---|---|
| `the_two_tables_agree_cell_for_cell` | the windowed table folded over its windows equals the read-group table, cell for cell, at one read group (spec §12.6) |
| `one_walk_and_four_shards_give_identical_tables` | the same territory as the catalog's regions, and as pieces dealt to four accumulators and merged, gives identical tables |
| `no_locus_overlaps_the_one_before_it` | `loci_overlapping_previous` is zero, and so is a second count made outside the accumulator |
| `the_generic_path_fits_a_real_sample_without_railing` | the coupled fit returns, converges, and lands on neither end of the error-rate ladder |

The fourth is beyond the plan's three identities and is **added deliberately**: Checkpoint F
reads *"step 4's generic path runs end to end on real alignments, **and** the three
structural identities hold"*, and nothing else in the milestone runs a fit over loci a walk
produced. `F` is supplied rather than fitted, because both cohorts' BEDs hold a few hundred
windows against `MIN_WINDOWS_TO_FIT_INBREEDING`'s 3,000.

## The runs — five alignments, twenty test instances, all green

`tmp/run_f3.sh`, one invocation per row:

| alignment | generic loci | positions | reads | occupied cells | sites at the cap |
|---|---|---|---|---|---|
| HG002 30x BAM | 551,844 | 552,284 | 16,618,807 | 181 | 0 |
| HG002 300x CRAM | 550,049 | 552,625 | 67,982,188 | 155 | **545,863** |
| tomato SRR7279481 | 7,424,484 | 7,429,336 | 77,080,043 | 322 | 0 |
| tomato SRR7279482 | 7,348,533 | 7,359,219 | 194,900,571 | 377 | 0 |
| tomato SRR7279483 | 7,213,401 | 7,224,396 | 103,069,962 | 565 | 9,273 |

`loci_overlapping_previous` is **zero on all five**, and so is the independent count. The
sharded arm cut 3,142 HG002 typed regions into 4,773 pieces and 41,824 tomato ones into
63,357, over four accumulators, with every table identical.

## The one run that matters most, and it was nearly not made

**The 300x HG002 CRAM is the only alignment in either cohort that reaches the depth cap** —
545,863 of its 550,049 sites are subsampled to 124 reads, against zero on the 30x arm and
zero on two of the three tomato samples. That matters because **above the cap the two tables
are filled by different draws**: `count_by_read_group` subsamples each group against its own
depth, `count_whole_site_by_library` makes one shared draw for the site
(C3's amendment to C2). Their coinciding at a single read group is exactly what identity 1
claims, and no shallower run tests it. A first pass ran only the 30x arm, where the cap fires
nowhere.

The cap is still shipping without a *bias* measurement — no harness world reaches it (arch §8's
remaining `OPEN:`). What F3 adds is that it fires on real data, at 99.2% of sites on a 300x
sample, without breaking the identity the two tables rest on.

## Six mutations, six killed

Run one at a time against the HG002 30x arm, each reverted before the next; the tree was
byte-restored from a copy afterwards and `git diff` confirmed clean.

| mutation | outcome |
|---|---|
| `multi_library` forced **false** | killed by two *existing* accumulator tests (`a_multi_library_sample_keeps_the_attribution_in_the_windowed_table`, `merging_shards_that_disagree_about_the_library_count_is_refused`) |
| `multi_library` forced **true** | killed by identity 1 — the two tables then hold the same 181 cells with different keys |
| the read-group table charged `Bp(1)` per locus instead of the region's length | killed by identity 1's covered-position comparison |
| `merge` drops the read-group table | killed by identity 2 |
| `merge` drops the windowed table | killed by identity 2 |
| `in_pieces` returns its input | killed by identity 2's premise assertion |

**The `multi_library` question F1's review left open is now half-settled.** A reviewer's probe
had found two read groups over the same sites returning equal parameters whether the flag was
set or not, which raised the question of whether the flag is load-bearing at all. It is: three
tests die when it is flipped in either direction. What is *not* settled, and cannot be from
these cohorts, is whether the attributed key changes the **fitted numbers** on a genuinely
multi-library sample — every sample in both cohorts carries one read group, so F3 never enters
the attributed arm at all. That needs a multi-library alignment, and neither cohort holds one.

## Two defects in this work, both mine, both found by running it

**The rail check was inverted and called a perfectly ordinary fit railed.** The error-rate
ladder ascends in *Phred* and therefore descends in error rate — rung 0 is Phred 10, or 0.1,
and the last rung is Phred 50, or 10⁻⁵ — so `first()` is the coarsest rate, not the lowest.
Naming the ends `lowest`/`highest` by position rather than by magnitude made
`value > lowest && value < highest` reject HG002's Phred-26.5 fit on the first run. The names
are now `coarsest`/`finest` and the ladder's direction is stated where the comparison is made.

**The sharded arm's piece width did nothing.** At a fixed 10,000 bases it left HG002's 3,142
typed regions as 3,142 pieces — region typing has already fragmented the BED's spans to a mean
of 176 bases, so not one reached the width — and the arm therefore compared a walk against the
*same* region boundaries and tested only the merge. It passed, and said so nowhere. `in_pieces`
now cuts every generic region into **at least two** pieces whatever its length, with the width
as a ceiling on top of that, and the test asserts `pieces > regions` before walking. Nothing in
the earlier version could have reported the no-op; the number was visible only because the
`eprintln!` happened to print both.

## What F3 does not establish

- **It is a plumbing check, not evidence about the four numbers.** Spec §1 says so of identity
  1 in as many words: both tables reduce the same locus through the same counting functions and
  bin it with the same shared ladder. Identities 2 and 3 compare a walk against itself, and
  identity 4 asserts only that the fit returned and did not rail. The values are Milestone G's.
- **A zero counter cannot see a silent detector.** `loci_overlapping_previous` is counted twice
  — once by the accumulator and once over the loci themselves — so the *partition* claim rests
  on two implementations. Neither can catch a detector that never fires, because on data with
  no overlaps a deleted comparison answers zero too. C3's overlapping-loci fixture is what pins
  the positives, and it is a fixture precisely because no real cohort supplies the condition.
- **The attributed arm of the cell key is never entered**, as above.
- **`F` is never fitted.** 3,000 windows is 300 Mb that must hold sites; the tomato BED is 8 Mb
  and the GIAB one 572 kb.
- **The sharded arm cannot see a locus the generator drops at every boundary alike.** It can see
  one generated differently at a seam, and any count `merge` fails to carry.

## An observation for Milestone G, asserted nowhere

The 30x BAM and the 300x CRAM are the same reads — the first is the second downsampled with
seed 42 — over the same 100 spans, so the pair is two rungs of G2's coverage sweep arriving
early:

| | error rate | ladder rung | heterozygosity | hom-non-ref |
|---|---|---|---|---|
| HG002 30x | 2.239 × 10⁻³ | 66 (Phred 26.50) | 0.001407 | 0.000533 |
| HG002 300x | 2.371 × 10⁻³ | 65 (Phred 26.25) | 0.001543 | 0.000526 |

The error rate moves **one rung** and the homozygous-non-reference rate 1.3%, both inside what
the design argues a caller cannot feel. **Heterozygosity moves 9.7%, upward with depth**, which
is the direction a shallow arm losing heterozygotes would produce and the slope G2 exists to
name. Two points are not a sweep and nothing here is asserted; it is recorded so G2 knows where
to look first.

## Validation

All via `./scripts/dev.sh`, at the F3 commit:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --all-features` — **3,209 passed, 0 failed, 9 ignored**
  (the four new tests are `#[ignore]`d, so the passed count is unchanged and the ignored
  count rises from 5).
- `cargo test --lib parameter_estimation -- --list | grep -c ': test$'` — the module holds
  **311** tests, four of them ignored. **F2's report and its handoff both say the module held
  308 before this step; it held 307** — counting `#[test]` over
  `src/ng/parameter_estimation/` at `5d7c9a6e` gives 307, and so does the `--list` command
  once F3's four are subtracted. One test out, and again a number about the author's own
  tests rather than one copied from a design document.
- `cargo doc --no-deps --lib` — 12 unresolved links, the pre-existing baseline, none in this
  module.

`cargo test --all-targets` is still red through `benches/psp_writer_perf.rs:386`, in frozen
`src/psp/` code, ruled out of scope by the owner.
