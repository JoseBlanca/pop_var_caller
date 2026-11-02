# The hidden-duplication filter at low depth: σ₀ missed the reads' own noise

*2026-11-01, branch `dup-filter-depth` from main 0e635332. Found while evaluating neighbour-evidence
models on tomato1 (`doc/devel/ng/research/neighbour_evidence_hmm.md`, on branch `neighbour-hmm`),
where the filter dropped 95 records in 100 at one read a position. The filter is
`doc/devel/ng/spec/hidden_paralog_filter.md`; the coverage model it scores depth with is §3.1 there.*

## 1. What was wrong

The filter asks, per record, whether each sample's depth around the site looks like one copy or
several. "Looks like" is judged against σ₀, the scatter of a one-copy window's relative depth,
fitted per sample from the depths of its 500-base windows. **σ₀ was too small at low depth, and the
reason is one sentence of arithmetic.** `window_coverage.md` §4 said a window of 500 positions at
three reads a position "is a mean over ~1,500 reads". It is a mean over about ten: 1,500 is the read
*bases*, and with 150-base reads each read covers up to 150 of the window's positions, so the
positions share their reads. The fit's own comment rejected a sampling floor on the same premise
("a window's depth is far less noisy than a single Poisson draw").

The fit then measures the spread only inside a band of 0.4 to 1.6 times one copy, and at low depth
the band cuts away the very tails the reads' noise produces. On one tomato accession
(SRR7279540) at three depths, the spread of its 500-base window depths, relative to the median
(samtools depth over tomato1's 80 regions):

| depth | spread over all windows (1.4826 × MAD) | inside the band | windows outside the band | pvc's fitted σ₀, cohort median |
|---|---|---|---|---|
| 0.9× | 0.78 | 0.42 | 44 in 100 | 0.22 (every accession at 1×) |
| 2.9× | 0.45 | 0.34 | 20 in 100 | 0.30 (every accession at 3×) |
| 31× | 0.19 | 0.17 | not counted | 0.23 (native depth, 2.6× to 30×) |

At 1× every sample's ordinary depth noise then read as extra copies. Over 57 accessions, 106,238 of
118,405 records scored more than 5 natural-log units towards a duplication, the run-wide duplication
rate fitted from those scores came out at 0.95, and the cut dropped 113,041 records — 95 in 100,
including records whose own evidence favoured a real variant 140 to 1.

## 2. The change

σ₀ is the band's spread plus, in quadrature, the variance a window's mean depth has from which
reads happened to land in it. For reads of span `L` starting independently in a window of `W` bases
at one-copy depth `d`, summing each read's squared overlap gives, relative to `d²`:

    L ≤ W:   L·(W − L/3) / (W²·d)        — 0.27/d for 150-base reads in a 500-base window
    L > W:   (L − W/3) / (L·d)            — long reads, where the mean is about one Poisson count

`L` is measured, not configured: the window accumulator sums each sample's reported depth and the
rises in that depth from one position to the next (a read adds one where it starts), and their ratio
is the mean read span (`CoverageByGcHistogram::mean_read_span`). Both modes feed the accumulator the
same depths, so both get the same span, and nothing new is stored in a psp. Where a read ends at the
step another starts the rise is lost, so at high depth the span comes out long; there the term is
a thousandth of the band's variance and the error cannot matter.

The term is 0.27 at 1×, 0.009 at 30× and 0.0009 at 300×, against about 0.17² = 0.03 for the band on
a deep sample: it acts where depth is low and fades with depth.

## 3. What it does

### 3.1 tomato1, 57 accessions

Reference: the records the full-depth run's filter (main, without the change) kept — the only
verdict available on this cohort, not an independent truth. "Conditions" cut accessions to lower
depth with `samtools view -s` and call from psps with each condition's own fitted parameters.

| condition | σ₀, main → change | duplication rate | records dropped | positions the deep filter kept, dropped |
|---|---|---|---|---|
| every accession at 1× | 0.22 → 0.55 | 0.95 → 0.15 | 95.5% → 12.5% | 95.8% → 2.0% |
| every accession at 3× | 0.30 → 0.43 | 0.17 → 0.12 | 15.5% → 11.6% | 4.9% → 1.0% |
| nine deep accessions at 3×, the rest native | 0.25 → 0.30 | 0.13 → 0.12 | 12.3% → 11.7% | 0.9% → 0.6% |
| native depth (2.6× to 30×, median about 11×) | 0.23 → 0.27 | 0.118 → 0.115 | 11.7% → 11.4% | 0.0% → 0.3% |

Main here is 0e635332, which already aims each window's rate at a typical window's rather than the
run-wide one (`ng_dup_filter_misses_2026-11-01.md` §7). That does not help at 1×: every window's
rate is wrong in the same direction, and the typical window's rate is itself 0.95.

### 3.2 HG002 chromosome 20, one sample

GIAB benchmark BAMs at 5×, 10× and 30×, `--defaults`, one thread, flagged records tagged and then
scored by vcfeval inside the confident regions with the filter applied and not:

| depth | σ₀, main → change | records flagged on the chromosome | true calls the filter removes | false calls the filter removes | F1, filter applied |
|---|---|---|---|---|---|
| 5× | 0.263 → 0.350 | 13,388 → 11,955 | 61 → 0 | 22 → 2 | 0.7990 → 0.7993 |
| 10× | 0.204 → 0.262 | 15,457 → 14,442 | 50 → 1 | 24 → 3 | 0.9488 → 0.9490 |
| 30× | 0.143 → 0.172 | 18,149 → 17,918 | 10 → 2 | 85 → 84 | 0.9925 → 0.9926 |

At 30× the filter keeps all but one of the false calls it removed and stops removing 8 of the 10 true
ones. At 5× and 10× it nearly stops acting inside the confident regions: about 20 false calls a depth
come back, against 50 to 61 true calls kept. This is the depth tail the other report found at 5×
(`ng_dup_filter_misses_2026-11-01.md` §3.4: 77 of 79 true SNPs removed on clean ground sat at 1.8
copies or more by chance) and left for a design change to the depth term (§5.3 there); this is one
such change, the "spread fitted honestly" family rather than a count model.

### 3.3 What is not measured

**How many real duplications the change lets through.** The one set that counts duplications caught
— the simulated genome with planted duplications, where the filter catches 764 of 1,401 artifacts
at 3× (spec §3.3) — and the GIAB trio at 300× and the 20-sample simulated cohort are not on this
laptop. The change flags 11% fewer records on chromosome 20 at 5× and 1% fewer at 30×, nearly all
outside the confident regions, where nothing here says which were duplications. A wider σ₀ makes
the depth term weaker evidence in both directions: on the cross-platform fixtures the records'
ratios moved from −16 to −23 towards −5 to −6, every record still far from a cut.

## 4. Tests and checksums

- New: the accumulator counts the read span over overlapping reads, a gap, a position reported twice
  and a contig change; the sampling variance's two branches meet; σ₀ adds it to the band's spread.
- The histogram line the mode-equivalence oracle compares gains the two counts.
- Three calls checksums in `cli/cross_platform_digests.rs` re-recorded: only `PARALOG_LR` moved, and
  in the 21-sample fixture the header's `lr_cut` (32.25 → 31.75); every record still passes with
  posterior 0; leaving the term out reproduces each old checksum exactly. The fitted-parameters
  checksums did not move.
- 5,008 library tests pass.
