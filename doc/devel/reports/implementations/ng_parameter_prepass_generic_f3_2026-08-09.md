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

## The runs — five alignments, twenty test instances, nineteen green

Each row is one invocation of the module doc's command with the reads and the BED swapped
(HG002 30x BAM, HG002 300x CRAM, and three tomato CRAMs, against the GIAB and tomato1 BEDs):

| alignment | generic loci | positions | reads | occupied cells | sites at the cap |
|---|---|---|---|---|---|
| HG002 30x BAM | 551,844 | 552,284 | 16,618,807 | 181 | 0 |
| HG002 300x CRAM | 550,049 | 552,625 | 67,982,188 | 155 | **545,863** |
| tomato SRR7279481 | 7,424,484 | 7,429,336 | 77,080,043 | 322 | 0 |
| tomato SRR7279482 | 7,348,533 | 7,359,219 | 194,900,571 | 377 | 0 |
| tomato SRR7279483 | 7,213,401 | 7,224,396 | 103,069,962 | 565 | 9,273 |

`loci_overlapping_previous` is **zero on all five**, and so is the independent count. The
sharded arm cuts 3,142 HG002 typed regions into 7,453 pieces and 41,824 tomato ones into
97,732, over four accumulators.

**Identity 2 fails on tomato SRR7279481, and the failure is a real defect in locus
generation** — the next section. The other nineteen instances are green.

## The failure: cutting a region loses a deletion's tail

On `tomato1/crams/SRR7279481.p1.bench.cram` the cut walk yields **7,424,467 generic loci
against the uncut walk's 7,424,484** — seventeen fewer, over seventeen fewer positions, none
gained.

The seventeen are one contiguous run, `SL4.0ch01:32,931,592–32,931,608`, and they open a
piece that starts at 32,931,592. The read `SRR7279481.13095921` is aligned at 32,931,402 with
CIGAR `116M91D22M5S`, so its **91-base deletion spans 32,931,518 to 32,931,608** and its match
resumes at 32,931,609 — the first position that is *not* lost. The piece boundary falls 74
bases inside that deletion. **A generic locus widened to a deletion's reference span is clipped
at the region's end, and its tail is picked up by no later region.** The positions are covered
at depth 4 to 6, so this is not a coverage gap.

Three things narrow it:

- **It is the cutting, not the sharding.** One stream over the same pieces into one
  accumulator loses the same seventeen, so `merge` and the per-shard readers are not involved.
- **It is boundary-density dependent.** At an earlier, coarser cut — 63,357 tomato pieces
  rather than 97,732 — this sample passed. The other four alignments pass at the finer cut.
- **The case is known to the locus-generation plan and its fixture does not reach it.**
  `locus_generation_pileup_generator.md`'s fixture 6 is *"a deletion across a region boundary
  — the halo, checked against the same fixture walked as one region."*

**Left failing on purpose.** Loosening the assertion, or coarsening the cut until it passes,
would bury a defect that changes which sites a sharded run analyses. Repairing it belongs to
locus generation, which is another plan's territory — the same rule the plan states for
`loci_overlapping_previous`. **This is the stop-and-ask Milestone F ends on.**

## The run where the cap fires, and what it is actually worth

**545,863 of the 300x HG002 CRAM's 550,049 sites are subsampled to the ladder's cap of 124
reads.** One tomato sample reaches the cap too, at 9,273 of 7,213,401 sites — about one site
in 780 — and the 30x arm and the other two tomatoes never reach it. So the cap fires on two of
the five alignments, not one.

**What that buys belongs to identity 2, not identity 1**, and the first version of this report
said the opposite. It claimed that above the cap the two tables are filled by *different
draws* — `count_by_read_group` per group, `count_whole_site_by_library` once for the site —
so that their coinciding was a two-implementation agreement. `count_whole_site_by_library` is
reached **only when the accumulator's `multi_library` is set**, and this same report says four
paragraphs later that F3 never enters the attributed arm. At one read group both tables end in
`CountedSite::capped(depth, alt, cap, seed_at(region))` with the same four arguments and the
draw is pure in them: one draw, evaluated twice, above the cap exactly as below it. Worse, the
claim that no shallower run tests it is false —
`accumulators.rs`'s `the_two_tables_agree_cell_for_cell_above_the_cap_too` drives 40 synthetic
sites at depths 300 to 963, asserts all 40 were capped, and compares the two tables cell for
cell. Its doc had the mechanism right; I wrote a different one into four documents.

**The property the deep arm really does carry is identity 2's.** `seed_at` draws from the
contig and the start position and nothing else, so a capped site must keep the *same reads*
however the genome is cut into regions. On the 300x arm 545,863 sites are drawn that way, and
identity 2 is the only test in this module that can see that fail.

The cap is still shipping without a *bias* measurement — no harness world reaches it (arch §8's
remaining `OPEN:`).

## Seven mutations, seven killed

Run one at a time against the HG002 30x arm, each reverted before the next; the tree was
byte-restored from a copy afterwards and `git diff` confirmed clean. The last was added in
response to the review and is the one that mattered — see below.

| mutation | outcome |
|---|---|
| `multi_library` forced **false** | killed by two *existing* accumulator tests (`a_multi_library_sample_keeps_the_attribution_in_the_windowed_table`, `merging_shards_that_disagree_about_the_library_count_is_refused`) |
| `multi_library` forced **true** | killed by identity 1 — the two tables then hold the same 181 cells with different keys |
| the read-group table charged `Bp(1)` per locus instead of the region's length | killed by identity 1's covered-position comparison |
| `merge` drops the read-group table | killed by identity 2 |
| `merge` drops the windowed table | killed by identity 2 |
| `cut_into_pieces` returns its input | killed by identity 2's premise assertion |
| `add_locus` returns immediately | **added after review** — killed by all four, identity 2 at its new emptiness premise |

**The `multi_library` question F1's review left open is now half-settled.** A reviewer's probe
had found two read groups over the same sites returning equal parameters whether the flag was
set or not, which raised the question of whether the flag is load-bearing at all. It is: three
tests die when it is flipped in either direction. **One** of those kills belongs to tests that
already existed — forcing the flag false is caught by two of the accumulator's own — and the
other two are identity 1, which this step adds. An earlier version of this report and of
`2d8c291e`'s message said two of the kills were pre-existing. What is *not* settled, and cannot be from
these cohorts, is whether the attributed key changes the **fitted numbers** on a genuinely
multi-library sample — every sample in both cohorts carries one read group, so F3 never enters
the attributed arm at all. That needs a multi-library alignment, and neither cohort holds one.

## What the review changed

Three agents on `2d8c291e`'s diff, in isolated worktrees. **The finding that matters is a
Blocker one of them demonstrated by mutation: identity 2 passed over an empty walk.** It
guarded its *regions* — pieces outnumber regions, more than one shard — and never its *loci*,
so two empty accumulators satisfied every comparison in it: two empty ploidy sets, two vacuous
per-ploidy loops, two empty key lists, two default counter sets. `WalkTally` did not rescue it,
because the tally is built from `locus.region` outside the accumulator and stays non-zero when
nothing is entered. Gutting `add_locus` failed identities 1 and 3 and panicked the end-to-end
test, and left this one green. It now asserts `walked_once.loci > 0` and
`whole.covered_positions() > 0` first, and the same mutation now fails all four. **The same
green comes without any mutation, from a BED whose spans hold no reads** — which is what
someone gets on their first attempt at a new cohort.

Two more, both real:

- **The end-to-end test asserted a property of a fit without asserting a fit happened.**
  `resolve_error_rates` hands back `DEFAULT_ERROR_RATE` = 0.001 with `Provenance::Defaulted`
  for a read group below `MIN_SITES_TO_FIT` with no sibling to lend, and 0.001 sits
  comfortably inside the ladder — so a sample whose every group defaulted passed the rail
  check while nothing was fitted. Unreachable on two cohorts of single-library samples with
  half a million sites each, and this file takes its alignment path from the environment. It
  now asserts `Provenance::FittedHere` per group.
- **`cut_into_pieces` and `deal_into_shards` are pure functions whose every caller is
  `#[ignore]`d**, so `cargo test` compiled them and ran neither. They now have five unit tests
  that need no data — and the first run of those tests **caught an off-by-one in the thirds
  rule I had just written**: at `ceil(L/3)` a four-base region gets a width of two and comes
  back in two pieces, not three, breaking the guarantee the function exists for at the one
  length nobody checks by hand. `L / 3` fixes it.

Also applied: `argmax_at_ladder_end` is dropped between the fit and `GenericSampleParameters`,
so no caller can read the bit arch §9 calls one of the two ways this estimator returns a
confident wrong number — **recorded for the owner rather than fixed**, because carrying it
changes a public type and reaches past this step. An unreachable empty-BED assertion removed
(`RegionSet` already returns `BedError::NoRegions`). The `.fai` verification joined *before*
the assertions rather than after, since `VerificationHandle`'s `Drop` is silent while
unwinding and a stale index is the first hypothesis a reader of a red identity needs ruled
out. Cell comparisons now name the first differing cell instead of printing both tables, which
run 155 to 565 cells. §12.6 also checked on the collapsed table a supplied `F` builds, and
arch §9's other half — no cell scored below its own alternative count — asserted on a table a
real walk filled. Renames: `Target` → `WalkInputs`, `where_` → `run_label`, `PIECE_BP` →
`MAX_PIECE_BP`, `SHARDS` → `SHARD_COUNT`, `in_pieces` → `cut_into_pieces`, `in_shards` →
`deal_into_shards`, `required` → `required_env_var`, `describe` → `one_line_summary`;
`cells_of` deleted.

**One reviewer finding rejected, with a measurement.** An agent asked me to restore the
`debug_assert!` justification for `--release`. I had already replaced it because I measured
it: the debug build **passes** this walk — it runs no production code, which is what trips the
assertion in `parity.rs` — and takes 128.7 s against release's 12.4 s.

## Three defects in this work, all mine, all found by running it

**The `--release` justification was copied from `parity.rs` and does not hold here.** That
file needs release because real paired-end data trips a reachable `debug_assert!` in
*production's* walker, which it runs beside ng's. Nothing in F3 runs production, and the debug
build completes: `no_locus_overlaps_the_one_before_it` on the HG002 30x arm passes in debug.
What release actually buys is speed — **12.4 s against 128.7 s**, 10.4 times — and the doc now
says that instead. Stating a failure mode a path does not have would have sent the next reader
hunting a symptom that does not occur.

**The rail check was inverted and called a perfectly ordinary fit railed.** The error-rate
ladder ascends in *Phred* and therefore descends in error rate — rung 0 is Phred 10, or 0.1,
and the last rung is Phred 50, or 10⁻⁵ — so `first()` is the coarsest rate, not the lowest.
Naming the ends `lowest`/`highest` by position rather than by magnitude made
`value > lowest && value < highest` reject HG002's Phred-26.5 fit on the first run. The names
are now `coarsest`/`finest` and the ladder's direction is stated where the comparison is made.

**The sharded arm's piece width did nothing.** At a fixed 10,000 bases it left HG002's 3,142
typed regions as 3,142 pieces. Region typing cuts a BED span into runs of generic sequence
between repeat tracts, and those runs are short — 3,142 regions carry 552,284 covered
positions between them — so none came near 10,000, and the unchanged piece count is the
evidence of that. The arm therefore compared a walk against the *same* region boundaries,
tested only the merge, and said so nowhere. Nothing in that version could have reported the
no-op; the number was visible only because the `eprintln!` happened to print both.

The first repair halved every region, and **halving was still not enough**: it gives every
region an interior *seam* but never an interior *piece*, one with a cut at each end, which is
the case where a generator's leading and trailing context both fall on a boundary.
`cut_into_pieces` now takes **thirds**, so every region of three bases or more yields an
interior piece whatever the organism, with `MAX_PIECE_BP` a ceiling on top; the test asserts
`pieces > regions` and that some region produced three pieces, before walking. That change is
what exposed the locus-generation defect above — at halving, tomato SRR7279481 passed.

## What F3 does not establish

- **It is a plumbing check, not evidence about the four parameters this step fits.** Spec §1
  says so of identity 1 in as many words: both tables reduce the same locus through the same
  counting functions and bin it with the same shared ladder. Identities 2 and 3 compare a walk
  against itself, and the end-to-end test asserts only that every rate was fitted from the
  sample's own sites and none sits on an end rung. The values are Milestone G's.
- **A zero counter cannot see a silent detector.** `loci_overlapping_previous` is counted twice
  — once by the accumulator and once over the loci themselves — so the *partition* claim rests
  on two implementations. Neither can catch a detector that never fires, because on data with
  no overlaps a deleted comparison answers zero too. C3's overlapping-loci fixture is what pins
  the positives, and it is a fixture precisely because no real cohort supplies the condition.
- **The attributed arm of the cell key is never entered**, as above.
- **`F` is never fitted, and it is scatter rather than size that prevents it.**
  `fit_inbreeding` tests `windows_holding_sites` against 3,000 — three thousand *separate*
  100 kb windows each holding a site. The tomato BED is 8.0 Mb in 80 spans and the GIAB one
  572 kb in 100, so between them they touch a couple of hundred windows.
- **Every site is scored at ploidy 2.** `ConstantPloidy` is the only `PloidyMap` this plan
  builds and nothing here overrides it. Both BEDs are autosomal, so the assumption is met — but
  a BED reaching a sex chromosome, or a polyploid organism, would have its cells labelled with
  genotype classes they do not have, and identities 1–3 would not notice, because they compare
  the same labels against each other.
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

The error rate moves **one rung**, which is the ladder's own quarter-Phred spacing and what the
design calls finer than a caller can feel. The homozygous-non-reference rate **falls** 1.3%, and
that one the design does not cover: `spec/parameter_prepass_generic.md` §4 records the adopted
twenty-bin depth ladder as costing 0.3% of that rate and a sixteen-bin ladder as **rejected** at
1.8%, so a 1.3% move between two depths of one sample sits inside the range the design treats as
worth choosing a ladder over. **Heterozygosity moves 9.7%, upward with depth**, which is the
direction a shallow arm losing heterozygotes would produce and the slope G2 exists to name. Two
points are not a sweep and nothing here is asserted; it is recorded so G2 knows where to look
first.

## Validation

All via `./scripts/dev.sh`, after the review fixes:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --all-features` — **3,214 passed, 0 failed, 9 ignored** in
  the library target. The five new unit tests over the region-cutting helpers raise the passed
  count from 3,209; the four real-alignment tests are `#[ignore]`d and raise the ignored count
  from 5. (The command's *total* across its eleven binaries is larger — 3,283 — and earlier
  reports on this plan have quoted the library line as though it were the total.)
- `cargo test --lib parameter_estimation -- --list | grep -c ': test$'` — the module holds
  **316** tests, four of them ignored. **F2's report and its handoff both say the module held
  308 before this step; it held 307** — counting `#[test]` over
  `src/ng/parameter_estimation/` at `5d7c9a6e` gives 307, and so does the `--list` command
  once F3's nine are subtracted. One test out, and again a number about the author's own tests
  rather than one copied from a design document.
- `cargo doc --no-deps --lib` — 12 unresolved links, the pre-existing baseline, none in this
  module.

`cargo test --all-targets` is still red through `benches/psp_writer_perf.rs:386`, in frozen
`src/psp/` code, ruled out of scope by the owner.
