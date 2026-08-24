# ng — the two averages of a read's error, measured on real reads

**2026-08-24**, branch `ng-prereq-closeout`. The measurement owed to the owner by
[`ng_calling_prerequisites_d2_2026-08-24.md`](ng_calling_prerequisites_d2_2026-08-24.md) §6 and by
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.2, which until today said the
size of the gap was "unmeasured".

---

## 1. The answer

**The two averages are 25 times apart on the tomato cohort and 44 times apart on the deep human
sample. They are not close, they will never be close, and the question is closed the other way from
the one the expectation pointed at.**

| | tomato, 63 accessions, 2.5× to 28.6× | HG002, 100 benchmark regions, 300× |
|---|---|---|
| read-positions measured | 5,485,730,235 | 172,616,054 |
| **geometric** mean of the per-read error | 5.982 × 10⁻⁴ (Phred 32.2) | 2.905 × 10⁻⁴ (Phred 35.4) |
| **arithmetic** mean of the same errors | 1.505 × 10⁻² (Phred 18.2) | 1.282 × 10⁻² (Phred 18.9) |
| arithmetic ÷ geometric | **25.2** | **44.1** |

**On the depth figure for tomato, because it is not the one this project usually quotes.**
These are the benchmark slices, measured over the 8,000,000 bases of
`benchmarks/tomato1/regions.bed`, and there they run **2.5× to 28.6×, median 10.7×** — deeper
than the cohort's genome-wide figure of about three reads a position. So read the tomato column
as *a 63-sample cohort at ordinary depth*, not as *the low end of the depth axis*. The low end
here is the 2.5× accession, and its ratio is 37.0 — the **highest** of the 63, so the gap does
not close as depth falls.

Per read group on tomato the ratio runs from **22.7** to **37.0**, median 24.4 over 63 accessions.
The lowest is `SRS3394549_SRR7279515`; the highest is `SRS3394559`. **No read group is anywhere
near one**, which is what "the two means agree" would have looked like.

**What that would do to a charged error.** The scale is
`fitted error rate ÷ mean minted error`, and it is added once per observation in log space — the
caller holds no per-read error probability to multiply, only `exp(q_sum / num_obs)` per allele per
library. Building the scale from the arithmetic mean instead of the geometric one divides every
charged error by 44.1 on HG002 and by 25.2 on tomato — **16.4 and 14.0 Phred**. The reads would be treated
as that much cleaner than the pre-pass measured them to be. So the correction of 2026-08-24 was not
a tidy-up; it is the difference between a scale that is right and one that is wrong by a factor of
forty.

## 2. Why they are so far apart, and it is not the chemistry

**Three quarters of the arithmetic mean is reads the walk deliberately silenced.**

A read's minted error is `max(ln ε_BQ, ln ε_MQ)` — the worse of what the instrument said about the
bases and what the aligner said about the placement. When two mates of a pair overlap a position,
the walk keeps one and silences the other by giving it base quality Phred 0
([`genome_walk.rs`](../../../../src/ng/locus_generation/pileup/genome_walk.rs), the mate-overlap
rule). Phred 0 is an error probability of **exactly one**, and such a read still counts as a read.

`ln 1 = 0`, so a silenced read adds **nothing** to the log sum and a **whole unit** to the
probability sum. Measured:

| | silenced read-positions | their share of Σ ε |
|---|---|---|
| HG002 at 300× | 1,616,935 in 172,616,054 — 9 in 1,000 | **73.1%** |
| tomato, pooled over 63 | 38,440,389 in 5,485,730,235 — 7 in 1,000 | **46.6%** (per accession: 31.6% to 52.1%) |

So on the human sample, **9 read-positions in 1,000 account for nearly three quarters of the
arithmetic mean**. The arithmetic mean is largely a measurement of how often a mate pair overlaps.
It is not a property of the library's chemistry, which is what the error-rate scale is for.

The geometric mean is not moved by them in the same way: a silenced read contributes zero to the
log sum, so it pulls the mean log error towards zero only through the count, not through the sum.

**This is a stronger argument for the geometric mean than the one the decision was taken on.** The
owner's reason was self-consistency — the model charges `exp(q_sum / num_obs)`, so the scale must
divide by that same quantity. That reason stands on its own. What the measurement adds is that the
alternative was not merely inconsistent, it was measuring something else entirely.

## 3. The identity the tool also checks, and what it found

The tool sums both shapes at the one place a read's own error still exists as one read's, then
folds every locus through the pre-pass's own accumulator
([`calibration.rs`](../../../../src/ng/parameter_estimation/generic/calibration.rs)) and prints
both read counts. **Two things broke that identity when the tool was first run, and both were
found by the check rather than reasoned about.**

**The walk has two paths that emit observations, not one.** The general fold through
`OpenPileupRecord::finalise` is the one the design documents describe; the ordinary-column fast
lane in [`fast_column.rs`](../../../../src/ng/locus_generation/pileup/fast_column.rs) builds
`SequenceObservation`s directly and never opens a record. A census hooked into the general fold
alone missed 4,165,737 read-positions of 172.6 million on HG002 — 2.4 in 100.

**The walk builds records in a halo beyond the region it was asked for**, which the generator then
discards by start position (`records_outside_region`). Counting at build time overcounted by
1,489,219 read-positions — 9 in 1,000. Totals are now held pending, keyed by the locus's start,
and move into the answer only when the generator keeps that locus.

**With both fixed, the two paths agree exactly**: 172,616,054 on HG002, and on all 63 tomato
accessions, with no locus left unruled on in any of the 64 runs. That is the site-set requirement
of §3.2 checked on 5.7 billion real read-positions rather than argued from the code.

## 4. Two corrections the measurement forces on documents

**The count is read-*positions*, not reads.** An observation contributes its `num_obs` at every
locus it appears at, so a 150-base read is 150 of these. Verified: 172,616,054 over 571,984 bases
is 301.8 a base, which is the sample's nominal 300×. `MintedReadErrors`'s doc had reasoned from
"a billion reads", which is low by about 150×, and concluded that an `i64` running sum had
four-hundred-fold headroom. Recomputed:

- a human genome at 30× is about 9.3 × 10¹⁰ read-positions, and at the mean log error measured
  here (8.145 nats) that is 7.9 × 10¹⁷ scaled units — an `i64` holds **twelve** such samples;
- **the same genome at 300× is 7.9 × 10¹⁸ on its own, 86% of `i64::MAX`.** One sample.

So the widening to `i128` that the previous step made on a four-hundred-sample argument was
correcting a defect an ordinary deep human run would have reached. The doc now says so.

**The per-position depth cap does divide the two site sets, by 2.7% at 300×.** §3.2 said it did
not, on the grounds that the histogram's thinning is quality-blind. That argument is right **per
site** and is not the whole of it: the cap removes reads only from deep positions, so it changes
how much weight each position carries. Measured on HG002 at 300×, where the fit sees 70,288,390 of
172,616,054 read-positions:

| | geometric mean of the minted error |
|---|---|
| as the accumulator counts it (nothing thinned) | 2.9055 × 10⁻⁴ |
| with each position thinned to the histogram's cap first | 2.9862 × 10⁻⁴ |
| ratio | **0.9730 — 2.7%, or 0.12 Phred** |

On tomato it is nothing at all: on the deepest of the 63 accessions (28.6×) 228,468,065 of
228,492,796 read-positions are under the cap — 1 in 9,200 above it — and the mean moves by a factor
of 1.0000. **This is
an owner's decision and is left open**, recorded in §3.2 with both options and their costs: thin
the accumulator at the same cap, or accept 3 parts in 100 at 300× on the argument that the
population the scale is applied to at calling time is every read. Nothing decides it until the
scale has a consumer.

## 5. What was built

**[`pileup/minted_error_census.rs`](../../../../src/ng/locus_generation/pileup/minted_error_census.rs)**,
new — measurement scaffolding in the idiom `column_census` already established in the same module:
a process-global table, off unless `PVC_MINTED_ERROR_CENSUS=1`, read by one example and by nothing
in the caller. It keeps Σ `ln ε`, Σ `ε`, the read count, and how many reads carried `ε = 1`.

Three call sites, all behind the env check hoisted out of their loops: the general fold's
`finalise`, the fast lane's grouping pass, and the generator's region clamp (which is what turns a
built locus into a kept one).

**[`examples/ng_minted_error_means.rs`](../../../../examples/ng_minted_error_means.rs)**, new — the
driver. One line per read group with both means, their ratio, both read counts and whether they
agree, the silenced-read attribution, and the depth-cap column.

```
PVC_MINTED_ERROR_CENSUS=1 cargo run --release --example ng_minted_error_means -- \
    <reference.fa> <regions.bed> <sample.bam|cram> [sample ...]
```

**One difference from the real pre-pass, stated because it changes the site set.** Every BED
interval is walked as a *generic* region: the typed-region stream reads a repeat catalog built
beside the reference, and neither benchmark reference has one (`$HOME/genomes` is mounted
read-only, so the tool cannot write one). Repeat tracts inside the intervals are therefore walked
through the generic generator here, where the real pre-pass routes them elsewhere. It widens the
site set; it does not change what is being compared, because both means run over the same reads
whichever reads those are.

## 6. Validation

All in the dev container.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo clippy --example ng_minted_error_means` | clean |
| `cargo test --lib ng::locus_generation::pileup::minted_error_census` | 4 passed, 0 failed |
| `cargo test --lib` | **4,187 passed, 0 failed, 14 ignored**, 925 s — `main`'s 4,181 plus the six tests added across both commits |
| `cargo doc --no-deps` | 23 unresolved links, 12 redundant-explicit-link-target warnings — `main`'s baseline, unchanged |

**The `cargo doc` number is itself a correction.** The previous step's report called the baseline
24 and "unchanged"; it is 23. The census module first added two of its own — a doc comment on the
`pub mod` line resolves in the *parent* module's scope, where the names it linked are not — and
those are gone.

## 7. What is still owed

- **The depth-cap decision** (§4) — the owner's, and it needs the read likelihood's consumer to
  exist before it can be answered on anything but principle.
- **A read group with a borrowed error rate.** A group standing on fewer than 10,000 sites gets
  the mean of the other groups' rates rather than its own, while its denominator stays its own
  reads. §3.2's sentence about one site set does not describe that case; the module doc now names
  it. A capture panel, or a minor library in a multi-library sample, reaches it.
- **Nothing beyond 300×.** The cap divergence grows with depth and 300× is the deepest sample in
  `benchmarks/`.
