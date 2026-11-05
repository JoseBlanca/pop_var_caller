# Performance Review: pipeline
**Date:** 2026-11-04
**Reviewer:** rust-performance-review skill (orchestrator, single reviewer, no fan-out)
**Scope:** pvc's whole pipeline on the tomato speed cohort; acted on the neighbour-evidence pass and the CRAM decode
**Verdict:** Apply the listed wins (applied on branch `perf-pass`); the fit needs an owner decision
**Hot-path evidence:** macOS `sample` and `samply` profiles, DHAT heap profile, `/usr/bin/time -l` on every run

---

## 1. Scope and constraints

**What was reviewed.** pvc's four steps on the cohort the benchmark suite times for speed: 62
tomato accessions at about three reads a position, over 80 regions of 100 kb (8 Mb), at 6
threads. Measured on main `88f93350` (the reads-in-draw merge), changes on branch `perf-pass`.

**Targets** (owner, 2026-11-04): wall time and peak memory; a large gain in one at a small cost in
the other is a win; a change that makes the code clearer is kept even for a small gain, and one that
complicates it for a small gain is rejected.

**Why no sub-agent fan-out.** The skill dispatches one measuring agent per category, each in its
own worktree. This laptop has six fast cores, other sessions kept its load between 2 and 11, and
the owner's rule is that measurements share at most six threads. One agent measuring in turn gave
comparable numbers; parallel agents would not have.

**Where the time and memory go, before any change** (main, load 5 to 7, so wall times run high):

| step | wall | peak memory | instructions |
|---|---|---|---|
| `generate-psps`, 62 processes 6 at a time | 86 s | — | — |
| `estimate-parameters` (the fit) | 690 s | 390 MB | 91,778 G |
| `call-from-psps` | 130 s: calling 38 s, the neighbour-evidence pass 90 s | 688 MB | 9,488 G |
| `call-from-alignments` | 126 s | 14.5 GB | 7,105 G |

Three things stood out: the fit is three quarters of the psp route; the pass is 70% of
`call-from-psps`; and `call-from-alignments` held 14.5 GB, against 9 GB in the suite's October run.

**Categories applied:** methodology (one change per measurement, builds alternated, load logged),
allocations and data layout (the CRAM reference buffers), hot loops and float portability (the
pass's walk and draw, the fit's integral).

## 2. Verdict

**Apply the listed wins** — five commits on `perf-pass`, each measured alone. Together:

- **`call-from-alignments` peak memory falls by a third** on the 62 samples (14.3 → 9.2 GB), and by
  about 130 MB a sample on a whole human chromosome 20, where it grew with the chromosome. Output
  byte-identical.
- **The neighbour-evidence pass is about a third faster**, and `call-from-psps` about a fifth:
  alternating the builds, the pass took 91 and 89 s on main against 64 and 56 s, the whole run 137
  and 128 s against 109 and 95 s; 13% fewer instructions. The output moves, within the noise of the
  pass's random draws (§5).

**The fit is left alone, and is the open question** (§3): it is 690 of about 900 seconds of the psp
route, its cost is spread over many small sites, and the one experiment run on it gained no wall time.

## 3. Measurement plan

What was run, so the numbers can be repeated (`tmp/perf/` in the `perf-pass` worktree):

- `pipeline.sh <label> [catalog psps fit call align]` — the four steps under `/usr/bin/time -l`;
  `PVC_SAMPLES=n` takes the first `n` samples, `PVC_BIN` another build.
- `ab.sh <binary> <label>` — `call-from-psps` on the 62 samples: wall, peak, instructions, the
  pass's own time from the log, and a digest of the VCF body for bit identity.
- `vcf_diff.py` — how many genotypes two runs write differently, and how many differ in alleles.
- `sim.sh <binary> <label> <cells>` — the read-phase report's simulated cohorts, scored against
  their truth for wrong links, phased heterozygotes and wrong genotypes.
- `profile.sh`, `line_profile.py`, `inline_lines.py` — `sample` for a function ranking; `samply`
  plus `atos -i` for a line ranking that sees through inlining. `dhat_peak.py` for a DHAT profile
  (`pvc` built with `--no-default-features --features dhat-heap`, commit `c738c14c`).

**For the fit, if the owner wants it pursued**, the profile says where to look and that no single
line is worth much (8 samples, `samply`, share of pvc's samples by innermost line of our code):

| where | share |
|---|---|
| `fit.rs:4532-4534`: each position's expected genotype counts added per sample, node and genotype | 13% |
| `information.rs` 623-735: the scores the curvature is built from | about 14% |
| `fit.rs:4600-4622`: the duplicated class's expected counts | about 7% |
| `fit.rs:4241-4243`: the product over samples at each quadrature node | about 6% |
| `libm` `log` | about 4% |

Swapping the node and sample loops at 4241 so the 16 nodes' products advance together (bit-identical)
cut instructions 2.8% on 8 samples and gained no wall time (84.4 s against 81.1 s, load 4); it was
rejected under the owner's rule and kept as `tmp/perf/rejected_fit_loop_swap.patch`. **A real gain
in the fit is probably not per-pass arithmetic** — two earlier reviews took those
(`perf_snp-fit_2026-10-18.md`, and the thinned starts of 2026-10-25) — **but fewer passes**: the
62-sample fit makes three starts of 36 to 54 passes. That is the step recorded as next on
2026-10-25 (a Newton step from the true curvature), stopped then by the owner.

## 4. Build / toolchain configuration

No change. `pvc` gained one opt-in build: `--no-default-features --features dhat-heap` makes DHAT
the allocator and writes `dhat-heap.json` when `main` returns (`c738c14c`). It is how §5's memory
finding was made, and it compiles out of every other build.

## 5. Code-level findings

### Hot-path

**H1: [container.rs](../../../../src/read/input/aligned_reads_reader/container.rs) — the CRAM decode's reference reader held its contig** — applied in `913f5b26`.

To decode a CRAM slice, the reader fetches the stretch of reference the slice covers and copies
it out. Its reference reader is a window that extends forward when the next stretch starts close
behind — which consecutive slices always do — and nothing on the CRAM path asked it to drop what
lay behind. So every open file's decode, and its prefetch thread, held its contig up to the walk's
position. **Found with DHAT**: on 8 samples over 10 regions, 860 of the 1,106 MB live at the peak
were these buffers. The decode now holds its reference as a `SliceReference` (raw bases and
eviction) and evicts before each slice.

| run | before | after |
|---|---|---|
| HG002 chr20, one sample, 5× CRAM | 340 MB | 207 MB |
| HG002 chr20, one sample, 30× CRAM | 373 MB | 247 MB |
| tomato, 8 samples, 80 regions | 2.20 GB | 1.60 GB |
| tomato, 62 samples, 80 regions | 14.3 GB | 9.2 GB |

Output byte-identical in all four; instructions unchanged. **Scale:** about 130 MB a sample on
chromosome 20, two copies of its 64 Mb; on chromosome 1 (248 Mb) that would be about 500 MB a
sample, 30 GB over 62 (inferred, not run). The test
`no_reference_a_cram_cursor_holds_grows_with_the_walk` fails without the eviction (one reader held
221,415 of a 400,000-base contig).

**What remains of the tomato run's 9.2 GB**, about 150 MB a sample, is not explained. The
benchmark's CRAMs were cut down to the 80 regions, so one slice holds reads megabases apart and its
reference stretch spans them; that is the likely cause, but it was not measured on whole-genome
tomato CRAMs, which this machine does not hold.

**H2: [model.rs](../../../../src/neighbour_evidence/model.rs) — the pass's walk: four changes** —
applied in `cbc562cf`, `71109de3`, `19dfc534`, `0c67c9a3`.

The pass walks each sample along each window five times, and each walk computes, at every site,
a probability for each of 289 template pairs (16 templates and a background, squared), forward and
backward, then draws nine paths from them. The profile put 78% of the pass's CPU in three
functions of that walk. Each change was measured alone on `call-from-psps`, 62 samples:

| change | instructions | the pass's time | output |
|---|---|---|---|
| `margins`: each row summed with one running total, the columns over the row's slice, so they vectorise | 9,491 → 9,211 G | its own profile share 26,249 → 4,595 samples | byte-identical |
| one `pair_transition` for the forward and backward steps, its coefficients taken once a site instead of dividing by `K` at every pair | 9,210 → 8,853 G | — | 3 of 10,173,208 genotypes differ, 1 in alleles |
| the paths share each site's totals (`SiteDraw`) | 8,853 → 8,230 G | 75 → 64 s | the paths change (below) |
| the pair sums in eight lanes (`lane_sum`, `lane_dot`), normalising by one inverse | level | 60, 60 → 55, 55 s, alternating | 5 genotypes differ, 2 in alleles |

**The draw.** Each path built all 306 weights of a site and scanned them; a round draws nine paths
from the same probabilities. Given the next state, a site's weights are a share of every state, of
the next pair's row and column, and of the next state itself, so the running totals are summed once a
site and each path picks a part, then a state inside it, in about 17 steps instead of 306. The draw
is exact — `a_sites_draw_follows_forward_times_transition` checks every state's share against the full
weights, and fails when one part is given the wrong share — but the uniform that picked a state now
picks another, so the paths differ. On the dense 64-founder simulation, over nine seeds each (a
temporary seed switch, not committed): wrong links 934 on average before (spread 890 to 1,028) and
961 after (862 to 1,048), 1.2 standard errors apart; wrong genotypes 33,750 and 33,764. Over eleven
simulated cohorts, wrong genotypes 61,183 before and 61,155 after.

**The lane sums.** A sum of 289 values added one after another waits for each addition before the
next; eight running totals added in a fixed order do not, and the order depends only on the length,
so the bits are the same on every machine and pool width (`float_portability.md`). The gain shows in
wall time, not in instructions.

**Bits.** Three of the four changes move output bits, the same on every platform: the posterior
and phase hashes of `refine/tests.rs` were re-recorded on macOS arm64 and checked on Linux arm64
in the dev container (§6). Linux x86-64 was not checked.

**Peak memory.** `call-from-psps` peaked at 657 to 664 MB on main and 674 to 708 MB on the branch
over two runs each, 2% to 8% more. The draw now holds all nine paths' states at once rather than
one path's; that it is the cause is a guess, not measured.

### Likely

None open.

### Note

- **The fit**: see §3. Not acted on.
- **`call-from-psps` peak memory** rose by 15 to 50 MB (§5, H2); small against the run's 0.7 GB.

## 6. Out-of-scope observations

- **The pass walks each site twice on a whole genome.** Its windows are a 2 Mb core with 1 Mb
  margins on each side, so every site lies in two windows. The benchmark's 100-kb regions are one
  window each, so none of the figures above shows it. The margin is a modelling choice; halving it
  would halve the pass's work wherever windows tile a contig.
- **`generate-psps`** took 86 s for 62 processes; not profiled.

## 7. What's already good

- **The pass's arithmetic is already shaped for a fixed order**: every sum the walk does is over a
  slice in index order, so changing one of them to lanes was a local edit with a known effect on
  bits, and the hash tests caught each move ([refine/tests.rs](../../../../src/neighbour_evidence/refine/tests.rs)).
- **The CRAM cursor already exposed `resident_bases`** for exactly this kind of check; the new test
  needed no new instrumentation ([ref_seq.rs](../../../../src/ref_seq.rs)).

### Author response convention

Address each finding by its identifier with one of: `applied in <commit>` / `experiment shows no
gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`.

---

## Addendum: a second round on the pass (branch `perf-pass2`, the same day)

A line-level profile of the merged build (`samply` and `atos -i`, `tmp/perf/pass_lines.sh`) put a
quarter of the pass in the posterior and an eighth in the emission terms, both running inner loops
over a site's alleles — two at a biallelic site — so every running sum waited on the one before.
Six commits, each measured alone on `call-from-psps`, 62 tomato samples, **all bit-identical** to
main:

| commit | change | instructions |
|---|---|---|
| `4159786d` | the allele loops over slices, no bounds checks | 8,227 → 7,665 G |
| `248d27ca` | each normalisation a loop over one slice, not a chained iterator | 7,665 → 7,317 G |
| `3d08c0b0` | inner loops across the 17 templates, not the alleles: shares laid out allele by allele, pair weights column by column | 7,319 → 6,594 G |
| `0484c48f` | the posterior's states read in row order (no wall-time change; kept as the clearer form) | 6,592 → 6,494 G |
| `474bc426` | the draw reads the forward pass's margins and sums running totals only when a path needs them | 6,494 → 6,321 G |
| `32738444` | the likelihood matrix built once a site, not once a path | 6,323 → 6,198 G |

**Together, alternating with main:** the pass 60 and 58 s against 40 and 39 s; the run 100 and 96 s
against 81 and 78 s; instructions −25%; peak memory level (654 to 683 MB). Since main before this
review, the pass has gone from about 90 s to about 40 s on this cohort.

**Rejected:** writing the posterior's division as a select rather than a branch — the instruction
count did not move, so the compiler already did it.

**What is left** is the arithmetic each step needs: the posterior's division by each state's
emission (about 14%), the transition (5%), the margins and sums, and each path's haplotypes
(5%, several allocations a site a path, a change of the output's shape to remove). The larger
levers are modelling choices — the window margins, which double the work on a whole genome, the
rounds and the templates — and need accuracy measured before any is moved.
