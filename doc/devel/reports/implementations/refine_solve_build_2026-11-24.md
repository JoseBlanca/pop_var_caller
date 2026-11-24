# Refinement's solve: four exact speed-ups built (2026-11-24)

Branch `refine-solve`, from main 16c6f67f, which already holds the window-memory cuts of
[refine_mem_build_2026-11-24.md](refine_mem_build_2026-11-24.md). It builds the four changes the
solve's speed review measured and kept, items 4.1 to 4.4 of
`doc/devel/reports/reviews/refine_solve_speed_opportunities_2026-11-24.md` — a review on branch
`solve-hunt` (commit b6de77b6), not on this one. The review's prototype was written before the
memory cuts; each change was written again here by hand on the new code, one commit each. None of
the review's unbuilt items (fewer template candidates, single precision, fewer extra paths, keeping
the forward probabilities, building each site's shares) is built.

**What a reader should take away.**

- On a drawn cohort of 200 samples at 3 reads a position, the run retires 29% fewer instructions
  (2.379 to 1.683 × 10¹²) and the stage takes 26% less time (26.2 to 19.4 s, four threads). At 800
  samples, where every window fits the walk parameters, 40% fewer instructions (5.710 to 3.427 ×
  10¹²) and 39% less time (61.9 to 37.9 s).
- Most of it is the first change, which only acts in windows that fit the walk parameters — the
  first windows of a run, until each sample has seen 5,000 informative records. Over a whole
  genome those are a few windows; the other three changes act in every window, and together save
  4 to 7% of the instructions (3.7% at 800 samples, 6.5% at 200).
- The output is byte for byte main's after every commit, on 13 drawn cohorts from 1 to 800 samples,
  1 to 30 reads a position, thinned margins and not, two to six alleles.
- Peak memory did not rise beyond the run-to-run spread: with mimalloc, the largest resident
  memory was 510 to 554 MB on main and 523 to 528 MB on the branch at 200 samples, 564 to 607 MB
  and 589 to 605 MB at 800. The bytes asked of the allocator (dhat) peaked at 275.7 MB on main
  and 275.3 MB on the branch, at 200 samples.

Measured on a laptop (macOS, arm64) other jobs shared; at most six threads.

## What changed

| Review item | Commit | What |
|---|---|---|
| — | 556ed0c5, 3008f6e7, (this commit) | `examples/solve_bench.rs`, the review's timing tool, ported from `solve-hunt`; mimalloc as `pvc` uses it, and dhat when built with `dhat-heap`. |
| 4.1 | 9a3847b5 | The parameter fit walks forward only. |
| 4.2 | c7f64fae | Emissions from definite alleles by an allele-pair table. |
| 4.3 | 6817735c | Genotype pairs found by allele count, not hashed. |
| 4.4 | c471edb6 | Template choice reads a sample's two copies in one pass; the read-link draw marks its kept-out sites once a round. |

### 4.1 — the parameter fit walks forward only

Until a sample has seen 5,000 of its informative records, each window fits the sample's two walk
lengths by walking it at each of twelve grid points and adding up each walk's log-likelihood
(pass spec §3.5). Each was a full walk — forward pass, backward pass with every block's forward
probabilities computed again, posteriors — of which only the log-likelihood was read, and that is
the forward pass's alone.

`model::walk_log_likelihood` runs the forward pass alone, and `SampleFit::fit_window` calls it. The
review's prototype copied the forward loop; here the loop is one function, `forward_pass`, which
the full walk (`walk_with_draw`) runs too, with a callback through which the walk keeps a block's
first forward probabilities. So the two cannot drift apart, and the number the grid reads has the
bits of `walk(...).log_likelihood`. The persistence between two sites is computed on the fly
rather than read from the walk's precomputed list; it is the same function of the same gap, so the
same bits.

**Memory:** none. `walk_log_likelihood` takes the thread's walk arrays, as a walk does, and uses
two slots of the forward probabilities and one site's emissions; it shrinks the arrays' lengths
but not their capacity, so the next walk on the thread allocates nothing new.

### 4.2 — emissions from definite alleles by a table

From the second round a walk weighs its states by definite alleles: each template's shares at a
site are one allele drawn for certain (phasing spec §2.3). Two templates carrying the same allele
then have the same blurred shares, so a pair state's emission depends only on its two alleles:
`A × A` numbers rather than `K²` (4 rather than 289 at two alleles and 17 templates).
`definite_emissions_into` computes those numbers with the general code's arithmetic in its order
and copies each to the states carrying its alleles; at two alleles the copy is a choice between a
row's two numbers, which the compiler vectorises. A site where some template's shares are not one
allele for certain — a template carrying nothing, or a share that is not 0 or 1 — takes the
general computation. Ported unchanged from the prototype, apart from the scratch fields' comments.

**Memory:** three scratch vectors in each walk's emission scratch, `A`, `A²` and `K` numbers: under
1 KB.

### 4.3 — genotype pairs by allele count

`GenotypePairs` kept each allele count's genotype pairs in a `HashMap`, and a walk asks for a
site's pairs several times a site. A vector indexed by the allele count replaces the map, an entry
empty until its count is first asked for. **Memory:** one empty vector (24 bytes) for each count
below the largest asked for.

### 4.4 — template choice, and the kept-out sites

- `choose_haplotypes` (rounds 2 to 5) added each other sample's two haplotypes as two passes over
  the compared records, each reading again whether the sample is informative there.
  `StretchSums::add_both_copies` adds both in one pass, each copy's sums in the order its own pass
  would add them. The haplotype reader stays a `dyn` closure: the prototype's switch to a generic
  one measured nothing (review §6) and is not ported.
- The read-link draw asked at every site of each path a round draws (nine in an averaged round)
  whether the site keeps the sample's reads out — a hidden or read-check record — by two binary
  searches. The sites are now marked once a round. Unlike the prototype, the marks are made only
  when the draw reads reads (`draws_with_reads`); otherwise nothing is allocated.

**Memory:** one byte a walked site, for the round of one sample, on the thread solving it — 57 KB
at the 57,000-site window of the tomato replay.

## Verification

### Tests

After each commit, `cargo test --release --lib neighbour_evidence` and `--lib read_phase` passed:

| After | neighbour_evidence | read_phase |
|---|---|---|
| main | 143 | 33 |
| 4.1 | 144 | 33 |
| 4.2 | 145 | 33 |
| 4.3 | 146 | 33 |
| 4.4 | 147 | 33 |

They include the recorded digests of posteriors and phase, and the walk against the Python
reference to 10⁻⁹, all unchanged. New tests:

| Test | Checks |
|---|---|
| `model::tests::the_log_likelihood_alone_has_the_walks_bits` | `walk_log_likelihood` equals `walk(...).log_likelihood` bit for bit: on the four reference windows, and on 200 sites (more than three backward blocks) with template weights, a background landing share, three-allele sites and definite emission shares; and on one site. |
| `model::tests::definite_emissions_have_the_general_bits` | At 2, 3 and 5 alleles and 17 templates, the table's emissions equal, bit for bit, the general computation's on the same alleles given as template shares; the table takes a definite site and refuses one where a template carries no allele. |
| `model::tests::genotype_pairs_are_kept_by_allele_count` | Counts asked for out of order; a count not asked for is unknown. |
| `templates::tests::both_copies_in_one_pass_sum_as_two_passes` | `add_both_copies` gives the sums of two `add`s, bit for bit, over several stretches with records where a sample is not informative. |

At the end, `cargo test --release --no-fail-fast` over every target: 5,378 passed, 2 failed, 9
ignored (main after the memory cuts: 5,374 passed; the four more are this branch's). The two
failures are the tests that expect a panic only debug builds raise
(`alignment::left_align_structured::tests::an_offset_past_the_reference_panics`,
`genetics::tests::lgamma_panics_on_non_positive_in_debug`), which fail the same way on main under
`--release`.

`cargo clippy --release --all-targets`: no warning in any file this branch touches. One it raised
in `examples/solve_bench.rs` (a loop over indices) is fixed in 3008f6e7; the rest are in files the
branch does not change.

### Output, byte for byte

Every line the stage writes was hashed (SHA-256) on 13 drawn cohorts, with main's build and with
the build after each commit, at six threads. The first eight use `examples/refine_mem_probe.rs`
(multi-allelic: 2 to 6 alleles in tomato25's shares, 300-base fragments), with fixed windows and
an unbounded budget; the last five use `examples/solve_bench.rs` (biallelic, paired 350-base
fragments, windows capped so the first fits the walk parameters).

| Cohort | Samples | Reads a position | Windows | Thinned margins | Lines | Hash, every build |
|---|---|---|---|---|---|---|
| probe | 1 | 3 | 100 kb core, 100 kb margins, 300 kb | 5 a kb | 21,252 | f9296fc337c2b1a0 |
| probe | 2 | 1 | 100 / 100 kb, 300 kb | off | 29,551 | ad159dba3a4fdabf |
| probe | 50 | 3 | 100 / 100 kb, 300 kb | 5 a kb | 49,570 | 99bdfad4222b0ce5 |
| probe | 50 | 30 | 100 / 100 kb, 300 kb | off | 49,664 | 2fa83e2b8984e5e7 |
| probe | 50 | 30 | 100 / 100 kb, 300 kb | 5 a kb | 49,664 | 777d3290998882d6 |
| probe | 50 | 1 | 150 / 150 kb, 600 kb at 20 a kb | 5 a kb | 11,866 | 33533a27f17b8438 |
| probe | 200 | 3 | 30 / 30 kb, 100 kb | 5 a kb | 16,586 | 64aa791407ddcb59 |
| probe | 800 | 3 | 20 / 20 kb, 60 kb | 5 a kb | 9,953 | c772af985d0fd910 |
| bench | 200 | 3 | at most 50 kb, 144 kb | 5 a kb | 11,328 | a92a7cccd5f36447 |
| bench | 800 | 3 | at most 30 kb, 48 kb | 5 a kb | 3,961 | 8abb9c5242a8d373 |
| bench | 100 | 1 | at most 50 kb, 144 kb | 5 a kb | 11,116 | 92a886ef10c08d8e |
| bench | 50 | 30 | at most 25 kb, 72 kb | 5 a kb | 5,517 | cbf9092d65b00218 |
| bench | 50 | 30 | at most 25 kb, 72 kb | off | 5,517 | cbf9092d65b00218 |

The probe's records are 160 a kb unless the row says otherwise; the bench's are 12 bases apart.
The last two rows give one hash because no margin of so short a run lies 50 kb beyond a core,
where thinning starts; the probe's rows 4 and 5 are the ones that tell thinning on from off. The
pass does not run on one sample in a real run; the probe forces it on, so row 1 exercises the
code, not a configuration a user meets. Since the budget does not bind, the windows are the same
in every build, so a change of the output could come only from the solve.

## Measurements

`examples/solve_bench.rs`, the review's two configurations, four threads, one run at a time
(`/usr/bin/time -l`; instructions retired and cycles from it, the stage's time from the example).
These runs used the system allocator, as the review's did (the example got mimalloc afterwards,
for the memory runs below); main and branch alike.

- **M200**: 200 samples, 12,000 biallelic records 12 bases apart (144 kb), 3 reads a position,
  windows of at most 50 kb: one window that fits the walk parameters and two that do not.
- **M800**: 800 samples, 4,000 records (48 kb), 3 reads, windows of at most 30 kb: two windows,
  both fitting.

**Instructions retired, after each commit** (× 10⁹; main's three runs agree to 0.01%, the last
commit's to 0.003%):

| | main | 4.1 | 4.2 | 4.3 | 4.4 | all four |
|---|---|---|---|---|---|---|
| M200 | 2,379 | 1,839 | 1,750 | 1,714 | 1,683 | |
| saved, share of main | | 540 (22.7%) | 89 (3.7%) | 36 (1.5%) | 31 (1.3%) | 696 (29.2%) |
| M800 | 5,710 | 3,636 | 3,549 | 3,498 | 3,427 | |
| saved, share of main | | 2,073 (36.3%) | 88 (1.5%) | 51 (0.9%) | 70 (1.2%) | 2,282 (40.0%) |

These match the review's prototype on the same configurations (M200 2.372 to 1.683 × 10¹², M800
5.699 to 3.424 × 10¹²): the port lost nothing.

**Time, main against the last commit, three paired runs each** (the stage's seconds; one-minute
load 5 to 6 on 18 cores throughout):

| | main | branch | change |
|---|---|---|---|
| M200 | 23.6, 26.6, 26.2 | 19.4, 19.6, 19.2 | −26% (medians 26.2 and 19.4) |
| M800 | 61.9, 62.0, 61.4 | 38.0, 37.9, 37.7 | −39% (medians 61.9 and 37.9) |

Cycles moved with them: M200 3.65 to 2.68 × 10¹¹ (medians), M800 8.54 to 5.14 × 10¹¹. Single runs
after the middle commits: M200 18.9, 18.7 and 18.8 s; M800 41.0, 40.4 and 40.3 s.

**Where the first change acts.** Only in windows that fit the walk parameters. Both configurations
are short runs dominated by such windows — M800 entirely — so they show it at its largest. Over a
whole genome at 3 reads a position the fitting windows are the first one or two of each run (on
tomato, under 100 kb of records); more at 1 read, where about two records in three are informative,
and in a sparse genome. For a run limited to a region, a replay of a slice, or a test, it is the
largest part of the gain. The other three changes act in every window: on M200 they save 155 × 10⁹
instructions, 6.5% of main's run; on M800, 209 × 10⁹, 3.7%.

### Peak memory

The added memory is small by construction (each change's **Memory** line above: under 1 KB of
scratch a thread, and one byte a walked site for one sample's round). To check that nothing else
moved, `solve_bench` was built with mimalloc, as `pvc` is, for main's refinement sources and for
the branch's, and run three times each, interleaved, at four threads; the largest resident memory
from `/usr/bin/time -l`, in MB:

| | main | branch |
|---|---|---|
| M200 | 510, 554, 524 | 523, 527, 528 |
| M800 | 607, 564, 592 | 589, 593, 605 |

The two builds' ranges overlap, and each spans about 40 MB between runs of one build: no rise
shows. **dhat**, which counts the bytes the program asks the allocator for and so does not move
with the allocator's own choices, peaked at 275.7 MB on main and 275.3 MB on the branch (M200,
four threads, one run each; a run under dhat takes about 12 minutes, so M800 was not run).

The first timing runs, with the system allocator, also report macOS's "peak memory footprint",
which moved more between runs of one build than between builds (M200 main 488 to 615 MB, branch
486 to 509 MB); at M800 the branch's three runs (647 to 664 MB) sat above main's three (542 to 565
MB), while their largest resident memory did not differ (main 968 to 1,011 MB, branch 949 to 987
MB). The mimalloc and dhat runs above were made to settle that.

## Assumptions and choices

- **The forward loop is shared, not copied.** The prototype's `walk_log_likelihood` repeated the
  walk's forward loop; the review named the risk that the two drift apart. Here both run
  `forward_pass`, passing the site's persistence through a closure (a dynamic call a site, next
  to the `K²` work of the step) and a callback for the walk's checkpoints.
- **The kept-out marks only when the draw reads reads.** The prototype allocated them for every
  phased round; a round that does not read reads in its draw never asks.
- **`GenotypePairs::known` of a count whose pairs are empty** answers `None`, where the map
  answered an empty slice. Only zero alleles have no pairs, and no record has zero alleles.
- **Not ported:** the prototype's switch of the haplotype reader from a `dyn` closure to a generic
  one (review §6: no gain), and its `walk_bench` example.

## Not built, for the owner

The review's §5 lists what is left, ranked; none is built. The two largest change the output or
the memory and are the owner's to rule on: template choice in rounds 2 to 5 among fewer candidates
(most of a quarter of the solve at 800 samples, more at 2,176), and the walk in single precision.

## Reproducing

```
cargo build --release --example solve_bench --example refine_mem_probe
RAYON_NUM_THREADS=4 /usr/bin/time -l target/release/examples/solve_bench 200 12000 12 3 50 1 out.vcf
PROBE_OUT=out.vcf target/release/examples/refine_mem_probe 50 300 160 3 100 100 5
```
