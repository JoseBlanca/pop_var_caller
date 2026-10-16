# D3: the production repeat-tract fit at kimura's size, drawn (2026-10-16)

Plan step D3 of [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md) (spec §2.2, §4
item 5). Branch `str-penalised-fit`; the bench was built from commit `8c6ce6f5`, whose fitted numbers
are bit-identical to every later commit's fitting code up to this report.

## 1. The answer

**The production fit takes 1.4 to 2.9 minutes at 18 threads and 2.5 to 4.2 minutes at 6, with a peak
of 10,801 MiB (10.5 GiB), on a drawn cohort of kimura's size: 2,169 samples and 66.2 million (locus, sample) pairs
with reads.** The budget is one hour on one machine (research plan §6). The measurement is of the fit
alone; reading 2,169 censuses and the SNP/indel fit come on top and were not measured here.

| cohort variant | threads | fit, summed over the six motif lengths | largest process peak |
|---|---|---|---|
| `depth3`: 3 reads a sample, 32,118 loci | 18 | 83 s | 10,801 MiB |
| | 12 | 116 s | 10,610 MiB |
| | 6 | 151 s | 10,355 MiB |
| `depth1.35`: kimura's 41,257 loci, 74% of pairs with a read | 18 | 173 s | 10,654 MiB |
| | 12 | 193 s | 10,555 MiB |
| | 6 | 252 s | 10,427 MiB |

**Thread scaling.** From 6 to 18 threads the fit runs 1.8 times as fast on `depth3` and 1.5 times on
`depth1.35`. This Mac has 6 performance cores and 12 efficiency cores, so threads 7 to 18 run on
slower cores; the figures say how this machine scales, not how the code would on 18 equal cores. The
one-base motif length is 99.3% of the work; the small motif lengths (2 to 6 bases, under 0.4 million
pairs each) run slower at 18 threads than at 6 on `depth1.35` (two-base 14.5 s against 9.1 s): too
little work to fill the pool. They are seconds either way.

**What sets the time.** In the one-base process at 18 threads (`depth3`), 56.8 of the fit's 75.4 s
are the expectation step (the rounds of counting, 291 rounds), 11.2 s reading F off pairs of
samples, 6.8 s the standard errors.

**Memory.** The peak is the drawn evidence plus the fit's working memory; at 100 and 200 samples it
grew by 173 bytes a pair (`depth3`), of which the evidence is 155, and the full run's 10,801 MiB is
within 1% of that line's prediction (10,920 MiB). Free memory never fell below 34.9 GiB.

## 2. The cohort, and what was assumed

kimura's per-class counts were not available, so the per-class spread is assumed, built from the
tomato SL4.00 repeat catalog and kimura's totals (86,688 tracts in 141 classes, 41,256 kept by the
repeat-tract guard, 107 classes with fewer than 8 tracts with reads, 66.2 million pairs with reads).
The full derivation, its sources and the per-class table are in the preparation notes
(`tmp/d3/README.md`, `tmp/d3/models/*/n2169/classes.tsv`; scratch, not committed). In short:

- 141 classes over motif lengths 1 to 6, a cap of 5,000 tracts a class, every one-base tract kept by
  the guard and 0.57% of longer motifs' tracts, so that 34 classes keep 8 or more tracts and 107 fewer,
  as on kimura.
- The truth is the bench's tomato-like PCR stutter and allele-length curves.
- The bench gives every sample the same mean depth, so it cannot match kimura's loci, its share of
  pairs with a read and its reads a pair at once. **`depth3`** keeps the pair count and about three
  reads a pair (the time and memory drivers); **`depth1.35`** keeps kimura's loci and its 74% share of
  pairs with a read, at 1.8 reads a pair.

## 3. How it was run

- One bench process a motif length (`--period`), `--methods D_tract_prior` (the production setting),
  `--fixture-out` (fits every locus once, no held-out folds, no truth scoring), at 18, 12 and 6
  threads, one process at a time on an otherwise idle machine, each under the memory guard (free
  memory logged every 30 s, killed below 8 GiB; none was).
- `estimate-parameters` fits the motif lengths one after another in one process, so its tract-fit
  time is the sum above. It holds every motif length's evidence at once; the one-base process already
  holds 99.3% of it, so the bench understates the peak by under 0.1 GiB.
- Outputs: `tmp/d3/runs/20261016_031656_n2169/` (`summary.tsv`, `production_equivalent.tsv`).

## 4. What this does not cover

- The census reading and the SNP/indel fit. On 59 tomato2 samples `estimate-parameters` took 61
  minutes, 60 of them the SNP/indel fit, whose three starting points each stop at the 200-pass limit
  (`reports/reviews/perf_estimate-parameters_2026-10-15.md`, finding L5). At kimura's size that fit,
  not this one, decides whether the command meets its hour.
- kimura's real per-class spread and per-sample depth variation.
- Since the fit is far inside its budget, the two speed-ups that would break the arch's
  one-class-at-a-time memory rule (that review's L3 and L4) are not needed.
