# The duplicated class claims ten times more positions than the reads hold, and it picks them without looking at depth

*Research report, 2026-08-14, following
[`trio_heterozygosity_excess_2026-08-14.md`](trio_heterozygosity_excess_2026-08-14.md) §4.4, which
found the collapse on the tomato panel and said it needed its own investigation. Written for a reader
who has read none of the specifications.*

*Programs: `examples/ng_joint_records_walk.rs` — the 63 tomato accessions walked once, fitted with
the class off and with it on, and the claimed positions then measured against two things the fit
never looks at (`DUPLICATED_CLASS_AUDIT=<path>` turns the audit on). `examples/ng_joint_duplicated_drifting.rs`
— a drawn cohort whose allele frequencies come from drift rather than from the shape the fit assumes.
Raw output under `tmp/records/`, the scripts that produced it under `tmp/`.*

> **Note added when this report was brought onto `main` (2026-08-19).** Only the drawn-cohort
> program came across. **The audit of the 63 accessions cannot be re-run on `main` as written**: it
> reads each claimed position's depth against the sample's own GC-corrected median, and the
> coverage-by-window summary that supplied it was deleted from the estimator two days after this
> was measured (`refactor(ng): A3 — the coverage-by-window summary is deleted`). §4's depth test
> and §5's carrier counts therefore stand as a record of what was measured on 2026-08-14 rather
> than as something a reader can reproduce today. §6's drawn cohort is unaffected — that program is
> on `main` and runs. The per-position probability the audit reads,
> `JointFit::duplicated_posterior`, is on `main` too.

---

## 1. The question, and the answer

Before calling any variant, ng estimates how heterozygous each plant is. One thing corrupts that
estimate. Where a plant carries **two copies of a stretch of genome the reference holds once**, both
copies' reads pile up at the same position, and wherever the two copies differ from each other about
half the reads disagree with the reference — which is what a heterozygote looks like. The fit carries
a third class of position for exactly those, so they are not counted as heterozygous.

**Switched on, it takes the panel's median accession from 0.867 heterozygous positions per kilobase
to 0.064, and puts that accession 0.98 short of the heterozygosity random mating in the panel would
predict** — which says the plant is very nearly free of heterozygous positions. A rough SNP caller
measured the same cohort's pooled heterozygosity at 1.049 per kilobase. So the class claims about
**93% of the heterozygosity the pass exists to produce**, and nobody had measured which positions it
was claiming.

**It is over-claiming, and the reads say so two independent ways.**

- **Count.** The class books each accession as carrying an extra copy at **941 to 3,515 positions in
  every two million**, median 1,829. Counting the same accession's positions that actually look like
  a duplication — 35% to 65% of the reads disagreeing with the reference, inside a window carrying
  about twice that accession's normal depth — gives **24 to 468, median 183**, which agrees with the
  150 to 590 the duplication probe measured on eight alignments in
  [`duplicated_locus_probe_2026-08-12.md`](duplicated_locus_probe_2026-08-12.md) §6. **On the six
  accessions deep enough for the second count to be complete, the class books 8.9 times as many
  positions as the reads themselves single out.**
- **Depth.** The positions it claims sit at **ordinary read depth**. The median claimed position
  reads 0.86 to 1.36 times its accession's normal depth for that stretch's GC content, where a
  duplication reads 2.0. And the class does not use depth to choose: of the positions that read near
  half **in a window at twice normal depth** it takes 83% (median across accessions), and of those
  that read near half **at ordinary depth** it takes 86%. Those two numbers being the same is the
  whole finding.

**And a drawn cohort says how big the thing it is supposed to remove actually is.** Draw 63 samples
at tomato's depths and tomato's inbreeding, take the allele frequencies from **drift rather than from
the shape the fit assumes**, and plant duplications at the rate the probe measured on real
alignments. With the class off, heterozygosity comes back **11.3% above** what the cohort was drawn
with; with it on, 0.4% below. **The artefact is worth about a ninth of the heterozygosity, and on the
real panel the class removes thirteen-fourteenths of it** (§6.1).

**Of the three explanations, the evidence supports none of them cleanly, and what it rules out is the
first.** The class is not right (§3, §4). It is not specifically claiming rare variants either — its
claim rate *rises* with how many accessions carry the allele, from 1 position in 480 where no other
accession's reads show it to 1 in 2 where most of the panel's do (§5). And the population's frequency
density barely moves when the class is switched on, while the drawn cohort above shows the class
handling an out-of-family spectrum without difficulty, so it is not absorbing a misfit there either
(§6). **What it is doing is simpler than any of the three: it takes the heterozygous calls at every
depth and at every allele frequency**, because across a cohort a heterozygote and a duplication
carrier produce the same read counts, and the one thing that separates them — some accession being
homozygous for the non-reference allele — is a pattern this inbred panel has already made rare. Where
that pattern *is* present the class behaves: the most heterozygous accession in the panel keeps 74%
of its heterozygosity where the median one keeps 7.5% (§5.1). **The class is not broken; it is
starved, and on this cohort it is starved almost everywhere.**

**The obvious repair — hand the class the depth reading it has never had on this cohort — is not one
yet.** GC-corrected and fed to the fit, it brings the median accession back to 0.404 per kilobase,
about half of what the class-off fit books, and it does so by making the class claim **four times as
many positions**: 1 in 31 rather than 1 in 127 (§7).

**Recommendation: the class ships off.** §8 has the argument and what would earn the default back.

---

## 2. What was measured, and against what

**One walk, then fits that differ only in what the class is told, and everything else is arithmetic
on their output.** The 63 tomato bench accessions are walked once over
`benchmarks/tomato1/regions.bed` — 80 spans of 100 kb, of which 1,999,404 positions are kept — and
the resulting records are fitted with the class off, with it on (§3 to §6), and with it on and handed
each accession's own read depth (§7). The first two both converge in 30 passes, at 627 and 992
seconds.

The fit now keeps, for every kept position, **the probability that the position belongs to the
duplicated class**, and for every accession at every position **the probability that this accession
is the one carrying the extra copy** (`JointFitConfig::genotype_posteriors`, off by default at twelve
bytes a position an accession). Without those the class's claim is one number — a share — and a
share says how many positions it claims and nothing about which. With them the claim is a list of
positions, and a list can be checked.

**Two discriminators, and neither is a truth set.** There is no truth set for duplications in this
cohort, so everything below compares two ways of picking positions with each other:

- **Relative read depth**, which the fit does not use here. A window's mean depth is divided by what
  that accession's own coverage does at that stretch's GC content — on tomato the depth-against-GC
  curve spans a factor of 1.79, which is larger than the doubling being looked for — and rescaled so
  the accession's median window reads 1.0. Adjacent windows are summed until they hold about **12,000
  aligned bases**, which is 500 bp at 25 reads a position and 5 kb at 2.5, because below that a
  window's mean is scatter ([`duplicated_locus_probe_2026-08-12.md`](duplicated_locus_probe_2026-08-12.md)
  §4). Only windows that abut *in the genome* are summed.
- **How many accessions carry the alternative allele**, counted as *two or more of its reads show
  it*, which the fit uses only through the shape of the whole panel.

**The depth reading has a noise floor and it is visible in every table below.** Across the 63
accessions, **0.7% to 3.4% of all kept positions** read between 1.6 and 2.4 times normal — that is
what a position picked at random does, and it is what any claimed share has to beat.

---

## 3. The count, which is the cheapest check and already decides it

The class fits a share of **0.00789** of positions — 15,782 of 1,999,404, one position in 127 — with
the share of the panel carrying a given one drawn from `Beta(0.467, 2.701)`. **11,447 positions are
more likely in the class than not and 8,139 sit above a posterior of 0.9.**

Per accession, the class's own posteriors sum to **941 to 3,515 carrier positions per two million**.
Beside that stands the direct count, taken on the same kept positions in the same run: positions
where that accession reads 35% to 65% non-reference at eight reads or more, inside a window at 1.6 to
2.4 times its normal depth. **24 to 468 per two million.**

| | the class books, per 2 M | positions that look like a duplication, per 2 M | how many times more |
|---|---:|---:|---:|
| accessions at 20 reads a position or more (6) | 2,409 | 270 | **8.9×** |
| 10 to 20 reads (29) | 2,203 | 227 | 9.7× |
| 5 to 10 reads (18) | 2,044 | 175 | 11.7× |
| under 5 reads (10) | 1,541 | 60 | 25.5× |

**Read the first row and treat the last as an upper bound.** A position cannot show a near-half read
share until it carries several reads, so the direct count is deliberately restricted to positions
with eight reads or more — which at 2.4 reads a position is a small and unrepresentative slice of the
genome. The deep accessions are where both counts are over comparable sets, and there the class books
**about nine times** what the depth reading finds.

**The direct count reproduces the earlier measurement, which is what says it can be trusted here.**
Its median across the 63 accessions is 183 per two million against the 150 to 590 the duplication
probe measured on eight alignments over a different set of positions.

---

## 4. H1 — the depth test, per accession

**A stretch of genome an accession carries twice collects two copies' reads.** So if the class is
picking duplications, the positions it claims should read about twice that accession's normal depth.
Three columns say whether they do: the claimed positions, the positions the *same fit with the class
off* calls heterozygous — which is where the claims came from — and the whole genome.

**The median claimed position reads 1.08 times normal** (0.86 to 1.36 across accessions), where a
duplication reads 2.0. **Between 9.3% and 36.1% of claimed positions sit in a window at 1.6 to 2.4
times normal**, median 14.7%, against a background of 0.7% to 3.4%. So the claimed set is enriched in
two-copy windows, by roughly a factor of ten — and so is the heterozygous set it came from.

**In 53 of the 63 accessions the positions called heterozygous with the class off are *more* often in
a two-copy window than the positions the class claims**, median 18.0% against 14.7%. The class is
therefore not concentrating the duplications: the set it takes is very slightly *less*
duplication-like than the set it takes it from.

**The sharpest form of the same statement, and it needs no threshold on the class's side.** Take the
positions where an accession reads near half at eight reads or more, and split them by whether the
window is at two copies or at one:

| | across the 63 accessions |
|---|---|
| near-half positions in a window at 1.6–2.4× normal — the class takes | 63.4% to 95.1%, **median 83.4%** |
| near-half positions in a window at 0.6–1.4× normal — the class takes | 6.9% to 100%, **median 86.1%** |

**The class takes a near-half position at ordinary depth as readily as one at twice depth.** Whatever
it is selecting on, it is not the thing that separates a duplication from a heterozygote.

### 4.1 Every accession, deepest first

**The count, per accession.** `window` is how wide a stretch that accession's depth had to be read
over to collect 12,000 aligned bases. `looks duplicated` counts the positions where it reads 35–65%
non-reference at eight reads or more inside a window at 1.6 to 2.4 times normal.

| accession | reads a position | window | the class books, per 2 M | looks duplicated, per 2 M | how many times more |
|---|---:|---:|---:|---:|---:|
| SRS3394685 | 30.6 | 0.5 kb | 2839 | 228 | 12.5× |
| SRS3394606 | 27.0 | 0.5 kb | 2164 | 242 | 8.9× |
| SRS3394682 | 25.6 | 0.5 kb | 1625 | 61 | 26.6× |
| SRS3394674 | 25.1 | 0.5 kb | 2093 | 468 | 4.5× |
| SRS3394687_SRR7279542 | 24.4 | 0.5 kb | 2250 | 443 | 5.1× |
| SRS3394702 | 24.3 | 0.5 kb | 3485 | 176 | 19.8× |
| SRS3394714 | 19.3 | 1.0 kb | 2610 | 314 | 8.3× |
| SRS3394641 | 18.9 | 1.0 kb | 1762 | 171 | 10.3× |
| SRS3394701 | 17.9 | 1.0 kb | 3456 | 278 | 12.4× |
| SRS3394678 | 17.6 | 1.0 kb | 1627 | 151 | 10.8× |
| SRS3394696 | 17.2 | 1.0 kb | 2581 | 277 | 9.3× |
| SRS3394697 | 15.4 | 1.0 kb | 2622 | 244 | 10.7× |
| SRS3394707 | 15.3 | 1.0 kb | 1926 | 200 | 9.6× |
| SRS3394695 | 15.2 | 1.0 kb | 2662 | 287 | 9.3× |
| SRS3394604 | 14.3 | 1.0 kb | 2113 | 200 | 10.6× |
| SRS3394711 | 14.1 | 1.0 kb | 2013 | 205 | 9.8× |
| SRS3394692 | 13.8 | 1.0 kb | 1539 | 169 | 9.1× |
| SRS3394713 | 13.6 | 1.0 kb | 2496 | 222 | 11.2× |
| SRS3394605 | 13.1 | 1.0 kb | 1829 | 188 | 9.7× |
| SRS3394709 | 13.0 | 1.0 kb | 3432 | 309 | 11.1× |
| SRS3394699 | 12.6 | 1.0 kb | 1694 | 183 | 9.3× |
| SRS3394691 | 12.5 | 1.0 kb | 1792 | 193 | 9.3× |
| SRS3394710 | 12.3 | 1.0 kb | 2968 | 362 | 8.2× |
| SRS3394705 | 12.3 | 1.0 kb | 1750 | 202 | 8.7× |
| SRS3394703 | 12.0 | 1.5 kb | 1541 | 165 | 9.3× |
| SRS3394690 | 12.0 | 1.0 kb | 1743 | 194 | 9.0× |
| SRS3394626 | 11.7 | 1.5 kb | 2893 | 265 | 10.9× |
| SRS3394688_SRR7279539 | 11.4 | 1.5 kb | 1465 | 170 | 8.6× |
| SRS3394694 | 11.2 | 1.5 kb | 1639 | 139 | 11.8× |
| SRS3394549_SRR7279515 | 11.2 | 1.5 kb | 1537 | 179 | 8.6× |
| SRS3394700 | 11.1 | 1.5 kb | 2528 | 281 | 9.0× |
| SRS3394706 | 10.9 | 1.5 kb | 1774 | 185 | 9.6× |
| SRS3394693 | 10.6 | 1.5 kb | 1800 | 178 | 10.1× |
| SRS3394698 | 10.3 | 1.5 kb | 2590 | 352 | 7.4× |
| SRS3394599 | 10.0 | 1.5 kb | 3515 | 314 | 11.2× |
| SRS3394704 | 9.9 | 1.5 kb | 2115 | 357 | 5.9× |
| SRS3394689 | 9.8 | 1.5 kb | 3123 | 173 | 18.1× |
| SRS3394638 | 9.8 | 1.5 kb | 1715 | 222 | 7.7× |
| SRS3394712 | 9.5 | 1.5 kb | 1868 | 154 | 12.1× |
| SRS3394708 | 9.5 | 1.5 kb | 1568 | 160 | 9.8× |
| SRS3394610 | 9.3 | 1.5 kb | 2050 | 188 | 10.9× |
| SRS3394663 | 8.4 | 1.5 kb | 3112 | 101 | 30.8× |
| SRS3394686 | 8.2 | 1.5 kb | 2991 | 268 | 11.2× |
| SRS3394594 | 7.3 | 2.0 kb | 1533 | 127 | 12.1× |
| SRS3394598 | 6.7 | 2.0 kb | 1662 | 136 | 12.2× |
| SRS3394655 | 6.0 | 2.5 kb | 1673 | 272 | 6.2× |
| SRS3394632 | 5.6 | 2.5 kb | 1758 | 106 | 16.6× |
| SRS3394635 | 5.5 | 2.5 kb | 941 | 138 | 6.8× |
| SRS3394636 | 5.3 | 2.5 kb | 1480 | 234 | 6.3× |
| SRS3394689_SRR7279535 | 5.0 | 2.5 kb | 2816 | 110 | 25.6× |
| SRS3394642 | 5.0 | 2.5 kb | 1027 | 105 | 9.8× |
| SRS3394633 | 5.0 | 2.5 kb | 3255 | 102 | 31.9× |
| SRS3394611 | 5.0 | 2.5 kb | 2114 | 203 | 10.4× |
| SRS3394712_SRR7279484 | 4.9 | 2.5 kb | 1903 | 96 | 19.8× |
| SRS3394688 | 4.3 | 3.0 kb | 1803 | 80 | 22.5× |
| SRS3394560 | 4.3 | 3.0 kb | 1149 | 94 | 12.2× |
| SRS3394595 | 3.7 | 3.5 kb | 1057 | 75 | 14.1× |
| SRS3394549 | 3.4 | 4.0 kb | 1660 | 53 | 31.3× |
| SRS3394641_SRR7279529 | 3.0 | 4.0 kb | 1341 | 50 | 26.8× |
| SRS3394615 | 2.8 | 4.5 kb | 1440 | 24 | 60.0× |
| SRS3394687 | 2.7 | 4.5 kb | 1746 | 57 | 30.6× |
| SRS3394640 | 2.5 | 5.0 kb | 1350 | 33 | 40.9× |
| SRS3394559 | 2.4 | 5.0 kb | 1957 | 41 | 47.7× |

**The depth, per accession.** `median relative depth` is over the positions the class calls this
accession a carrier of at a posterior above a half, on a scale where 1.0 is its own normal for that
stretch's GC content. The three `in 1.6–2.4×` columns are the share of three different sets of
positions that sit in a window at two copies: the ones the class claims, the ones the class-off fit
calls heterozygous, and every kept position — the last being the reading's own noise floor. The final
two columns are the share of near-half positions the class takes, split by whether their window is at
two copies or at one.

| accession | reads a position | median relative depth | claimed, in 1.6–2.4× | heterozygous with the class off, in 1.6–2.4× | any position, in 1.6–2.4× | near half at 2×, taken | near half at 1×, taken |
|---|---:|---:|---:|---:|---:|---:|---:|
| SRS3394685 | 30.6 | 1.10 | 12.1% | 13.5% | 0.9% | 91.4% | 84.2% |
| SRS3394606 | 27.0 | 0.93 | 11.1% | 12.8% | 0.9% | 78.0% | 47.8% |
| SRS3394682 | 25.6 | 0.91 | 9.3% | 14.2% | 2.1% | 95.1% | 88.4% |
| SRS3394674 | 25.1 | 1.06 | 25.3% | 30.5% | 1.3% | 89.0% | 92.0% |
| SRS3394687_SRR7279542 | 24.4 | 1.05 | 22.6% | 27.2% | 1.3% | 88.3% | 88.2% |
| SRS3394702 | 24.3 | 1.16 | 13.7% | 14.7% | 2.4% | 92.7% | 80.9% |
| SRS3394714 | 19.3 | 1.13 | 13.9% | 17.0% | 1.4% | 79.7% | 88.7% |
| SRS3394641 | 18.9 | 0.97 | 12.3% | 8.8% | 0.9% | 88.6% | 23.2% |
| SRS3394701 | 17.9 | 1.14 | 15.7% | 15.7% | 2.1% | 75.5% | 50.4% |
| SRS3394678 | 17.6 | 0.93 | 10.8% | 14.5% | 0.9% | 84.6% | 86.9% |
| SRS3394696 | 17.2 | 1.12 | 15.1% | 18.8% | 1.4% | 78.1% | 82.1% |
| SRS3394697 | 15.4 | 1.09 | 14.0% | 4.9% | 1.9% | 63.4% | 6.9% |
| SRS3394707 | 15.3 | 1.04 | 13.7% | 9.5% | 0.9% | 79.5% | 20.1% |
| SRS3394695 | 15.2 | 1.20 | 20.0% | 12.1% | 1.5% | 90.5% | 20.6% |
| SRS3394604 | 14.3 | 0.93 | 12.0% | 17.9% | 2.2% | 73.7% | 90.1% |
| SRS3394711 | 14.1 | 0.99 | 12.3% | 17.1% | 0.9% | 86.3% | 88.5% |
| SRS3394692 | 13.8 | 0.99 | 14.9% | 21.0% | 0.9% | 84.9% | 87.4% |
| SRS3394713 | 13.6 | 1.10 | 13.3% | 15.3% | 1.6% | 90.6% | 83.5% |
| SRS3394605 | 13.1 | 0.98 | 12.2% | 17.0% | 0.8% | 81.9% | 84.3% |
| SRS3394709 | 13.0 | 1.20 | 16.0% | 14.3% | 1.8% | 84.1% | 82.4% |
| SRS3394699 | 12.6 | 0.86 | 14.3% | 13.1% | 2.8% | 76.5% | 31.3% |
| SRS3394691 | 12.5 | 1.01 | 12.0% | 17.3% | 0.8% | 78.4% | 84.7% |
| SRS3394710 | 12.3 | 1.24 | 15.4% | 18.0% | 1.4% | 79.7% | 81.8% |
| SRS3394705 | 12.3 | 1.04 | 13.0% | 17.8% | 0.8% | 78.5% | 86.8% |
| SRS3394703 | 12.0 | 1.00 | 13.5% | 20.1% | 0.7% | 75.7% | 87.9% |
| SRS3394690 | 12.0 | 1.00 | 12.1% | 17.6% | 0.8% | 78.9% | 89.4% |
| SRS3394626 | 11.7 | 1.15 | 13.9% | 15.7% | 1.5% | 85.9% | 84.2% |
| SRS3394688_SRR7279539 | 11.4 | 1.03 | 15.6% | 23.0% | 0.8% | 85.2% | 83.6% |
| SRS3394694 | 11.2 | 0.98 | 12.3% | 17.0% | 0.7% | 87.7% | 88.3% |
| SRS3394549_SRR7279515 | 11.2 | 0.86 | 15.9% | 22.5% | 1.7% | 79.8% | 92.5% |
| SRS3394700 | 11.1 | 1.15 | 15.4% | 19.4% | 1.4% | 83.4% | 89.8% |
| SRS3394706 | 10.9 | 1.01 | 11.2% | 16.6% | 0.7% | 77.2% | 86.1% |
| SRS3394693 | 10.6 | 1.01 | 10.2% | 16.3% | 0.7% | 77.0% | 89.5% |
| SRS3394698 | 10.3 | 1.12 | 25.7% | 24.7% | 1.6% | 71.2% | 41.2% |
| SRS3394599 | 10.0 | 1.18 | 15.8% | 15.3% | 1.9% | 86.9% | 83.5% |
| SRS3394704 | 9.9 | 1.05 | 21.7% | 28.2% | 1.3% | 86.2% | 87.6% |
| SRS3394689 | 9.8 | 1.12 | 11.8% | 13.7% | 1.4% | 82.8% | 84.6% |
| SRS3394638 | 9.8 | 1.33 | 17.6% | 20.7% | 0.9% | 81.8% | 84.4% |
| SRS3394712 | 9.5 | 0.94 | 15.4% | 19.7% | 3.4% | 91.9% | 86.6% |
| SRS3394708 | 9.5 | 0.89 | 15.0% | 21.4% | 1.6% | 83.1% | 85.2% |
| SRS3394610 | 9.3 | 1.06 | 15.1% | 19.9% | 1.0% | 87.1% | 80.6% |
| SRS3394663 | 8.4 | 1.13 | 10.5% | 11.7% | 1.4% | 94.8% | 86.2% |
| SRS3394686 | 8.2 | 1.22 | 16.3% | 19.3% | 1.6% | 87.7% | 84.5% |
| SRS3394594 | 7.3 | 0.97 | 14.7% | 24.0% | 0.8% | 73.1% | 85.2% |
| SRS3394598 | 6.7 | 1.01 | 14.1% | 15.7% | 0.8% | 76.9% | 45.0% |
| SRS3394655 | 6.0 | 1.25 | 29.5% | 35.9% | 1.3% | 85.8% | 87.9% |
| SRS3394632 | 5.6 | 1.01 | 11.3% | 17.7% | 0.8% | 74.6% | 87.5% |
| SRS3394635 | 5.5 | 1.24 | 32.8% | 47.5% | 0.9% | 78.7% | 87.3% |
| SRS3394636 | 5.3 | 1.36 | 36.1% | 43.2% | 1.3% | 85.7% | 94.4% |
| SRS3394689_SRR7279535 | 5.0 | 1.11 | 14.0% | 18.1% | 1.4% | 83.1% | 91.7% |
| SRS3394642 | 5.0 | 1.17 | 31.8% | 43.7% | 0.8% | 87.1% | 83.1% |
| SRS3394633 | 5.0 | 1.17 | 11.8% | 11.6% | 2.9% | 88.2% | 72.6% |
| SRS3394611 | 5.0 | 1.08 | 22.8% | 30.4% | 1.4% | 84.3% | 88.3% |
| SRS3394712_SRR7279484 | 4.9 | 1.01 | 14.4% | 18.2% | 3.0% | 89.5% | 88.5% |
| SRS3394688 | 4.3 | 1.05 | 12.9% | 18.5% | 0.7% | 84.8% | 77.0% |
| SRS3394560 | 4.3 | 1.15 | 27.7% | 40.5% | 0.9% | 82.9% | 89.4% |
| SRS3394595 | 3.7 | 1.21 | 33.5% | 48.5% | 0.8% | 90.7% | 100.0% |
| SRS3394549 | 3.4 | 0.98 | 13.3% | 20.4% | 1.2% | 83.0% | 85.5% |
| SRS3394641_SRR7279529 | 3.0 | 1.21 | 24.3% | 15.2% | 1.0% | 78.7% | 64.9% |
| SRS3394615 | 2.8 | 1.35 | 22.1% | 27.8% | 0.9% | 87.5% | 94.0% |
| SRS3394687 | 2.7 | 1.23 | 28.2% | 35.6% | 1.2% | 82.5% | 100.0% |
| SRS3394640 | 2.5 | 1.13 | 24.4% | 28.7% | 0.9% | 94.2% | 100.0% |
| SRS3394559 | 2.4 | 1.22 | 23.5% | 29.7% | 1.0% | 82.9% | 99.9% |

**Which accessions this test can arbitrate on.** All of them carry a depth reading, because the
window is widened until it holds enough aligned bases; what a shallow accession loses is not the
reading but its resolution — at 2.4 reads a position the window is 5 kb, so a duplication shorter
than that is diluted by the single-copy sequence around it, and the near-half count that the reading
is compared against is restricted to the few positions carrying eight reads. **The six accessions at
20 reads a position and above are where both halves of the comparison are complete**, and they say
the same thing as the rest: median relative depth 0.91 to 1.10, claimed share in a two-copy window
9.3% to 25.3%, and the heterozygous-with-the-class-off share the same or higher in all six.

---

## 5. H2 — how many accessions carry the allele

If the class were claiming rare real variants, its claim rate would climb as the allele got rarer.
**It does the opposite.**

| accessions carrying the alternative allele | positions | the class's mass there | share of them claimed | carriers the class sees | heterozygosity there, with the class off |
|---:|---:|---:|---:|---:|---:|
| 0 | 1,720,771 | 3,556 | 0.21% | 0.99 | 2,692 |
| 1 | 191,063 | 2,513 | 1.32% | 1.88 | 5,298 |
| 2 | 33,337 | 1,155 | 3.47% | 2.92 | 4,592 |
| 3 | 12,021 | 771 | 6.41% | 4.09 | 3,520 |
| 4–5 | 12,690 | 1,107 | 8.73% | 5.40 | 8,104 |
| 6–10 | 14,783 | 2,523 | 17.07% | 8.10 | 22,817 |
| 11–20 | 9,169 | 2,257 | 24.62% | 13.15 | 31,650 |
| 21–40 | 4,421 | 1,303 | 29.47% | 25.87 | 30,072 |
| 41–62 | 1,127 | 582 | 51.65% | 43.72 | 20,945 |
| 63 | 22 | 14 | 65.26% | 35.74 | 334 |

*An accession counts as carrying the allele where two or more of its reads show it. "Heterozygosity
there" is the number of accession-positions the class-off fit calls heterozygous inside that band —
130,024 over the whole census.*

**Three things to read off it.**

- **The claim rate rises with frequency**, from 1 position in 480 where no accession's reads show the
  allele, to 1 in 76 where one accession's do, to 1 in 2 where 41 to 62 of them do. That is the class
  working as designed: a real variant half the panel carries should leave some accessions homozygous
  for the non-reference allele, and a duplication leaves none, so the absence is most informative
  where the allele is commonest.
- **But nearly half the mass is where that argument has nothing to work with.** Of the class's 15,782
  expected positions, **7,224 — 46% — sit at positions two or fewer accessions carry**, and a variant
  one accession carries is expected to leave no non-reference homozygotes whether it is a duplication
  or not.
- **And the heterozygosity is where the claim rate is highest.** 105,818 of the 130,024 heterozygous
  accession-positions — 81% — sit at positions six or more accessions carry. **67.8% of all of it
  sits at the 11,447 positions the class claims outright.**

### 5.1 The class spares an accession in proportion to how heterozygous it already is

The same mechanism read from the accessions' side. Comparing each accession's heterozygosity with the
class off and on, **the median accession keeps 7.5% of it and 41 of the 63 keep under a tenth** — but
the keeping is not uniform:

| accession | heterozygosity, class off | class on | kept |
|---|---:|---:|---:|
| SRS3394697 | 4.648 per kilobase | 3.448 | **74.2%** |
| SRS3394641_SRR7279529 | 1.135 | 0.632 | 55.7% |
| SRS3394695 | 2.150 | 1.028 | 47.8% |
| SRS3394707 | 1.543 | 0.735 | 47.6% |
| … | | | |
| the median accession | 0.867 | 0.064 | **7.5%** |
| the accession that keeps least | 0.318 | 0.009 | 1.5% |

**The most heterozygous accession in the panel keeps three quarters of its heterozygosity and the
median one keeps a fourteenth.** That is the cohort discriminator working where it can: an outbred
accession's heterozygous positions sit at variants the rest of this inbred panel carries in the
homozygous state, and a position with non-reference homozygotes in it cannot be a duplication. Where
an accession's heterozygous positions are its own — rare, and carried by nobody else in a form the
panel can check — the same discriminator has nothing to look at, and the class takes them.

**A tomato panel makes this worse than a random-mating one would.** These plants are 0.79 less
heterozygous than random mating in the panel predicts, so a variant at a frequency of a half leaves
most of its carriers *homozygous* — and a position where many accessions all read about half
non-reference is genuinely odd for such a panel. The class is right that something is wrong at those
positions. What it has no way to decide is whether the something is a duplication in each of those
plants or a stretch of the reference that two parts of the genome both map onto, which is the
mismapped class's job and which the fit judges separately: it books 58,765 of the 1,999,404 positions
as more likely mismapped than not.

---

## 6. H3 — the frequency density, and a cohort drawn from outside the model

**Inside the fit, a population's allele frequencies are four numbers**: the share of positions
carrying only the reference base, the share carrying only a non-reference one, and the two shapes of
a Beta covering everything between. If a real population's frequencies are shaped differently from
what those four can express, a class free to claim positions at any frequency would soak up the
difference.

**On the real panel the density barely moves.**

| | positions carrying only the reference | the rest segregate with | the population's expected heterozygosity |
|---|---:|---|---:|
| the class off | 0.9629 | `Beta(0.555, 6.151)` | 4.895 per kilobase |
| the class on | 0.9702 | `Beta(0.564, 5.801)` | 4.153 per kilobase |

**The population's expected heterozygosity falls 15% while each accession's own falls thirteenfold.**
So the class is not rewriting what the population looks like; it is reassigning individual plants'
genotypes out of *heterozygous*, and the homozygote excess absorbs the difference by moving from
0.789 to 0.983. Nothing here is arithmetically inconsistent — it is the same data described as *a
diverse population of near-clonal plants* instead of *a diverse population of somewhat heterozygous
plants* — and only evidence from outside the model can say which description is right.

### 6.1 A cohort drawn from outside the fit's own family says the shape is not the problem

**A drawn cohort whose defects come from the model the fit assumes cannot say whether the model is
right**, which is why every earlier drawn control has been generous to the class. This one draws its
allele frequencies from **drift**: a segregating position's alternative allele is carried by `k` of
the panel's `2n` chromosomes with `k` drawn from `1/k` — the spectrum a neutral, constant-sized
population leaves in a sample, so most variants sit on one or two chromosomes — and those `k` copies
are then **dealt out among the panel's chromosomes**, so the accessions are not independent draws at
a shared frequency. The fit's Beta cannot express either property. Duplications drift the same way,
planted at **300 carrier positions per two million per accession**, which is the middle of the range
the duplication probe measured on real alignments.

Everything else matches tomato: 63 accessions from 2.4 to 30.6 reads a position, every plant 0.79
less heterozygous than random mating predicts, a read misreading at 0.0034 at an ordinary position
and 0.0254 at a mismapped one, and 1 position in 33 mismapped. 200,000 positions.

| the cohort, 63 samples | drawn het/kb | the class off | the class on | class weight it fits | drawn |
|---|---:|---:|---:|---:|---:|
| **drift, duplications at the probe's rate** | 1.223 | **+11.3%** | **−0.4%** | 0.00066 | 0.00050 |
| drift, no duplications at all | 1.218 | −0.2% | −0.7% | 0.00019 | 0 |
| the fit's own Beta, no duplications | 2.098 | +0.3% | +0.3% | 0.00000 | 0 |

*At 25 samples the same three rows read +7.4% / −0.8%, −0.9% / −0.9%, and −1.8% / −1.8%.*

**Two things follow, and together they close H3.**

- **The class works on drawn data even when the frequencies come from outside its family.** Handed a
  drift spectrum and duplications at the real rate, it takes heterozygosity from 11.3% above the truth
  to 0.4% below it, and the weight it fits is within a third of what was planted. So neither the
  spectrum's shape nor the class's arithmetic is what fails on tomato.
- **The artefact it exists to remove is worth about 11% of heterozygosity at the measured duplication
  rate.** On the real panel it removes 93%. **That is the size of the discrepancy in one comparison**:
  the same estimator, the same panel size, the same depths, the same duplication rate, and eight times
  more heterozygosity taken on real reads than on drawn ones.

**What the drawn cohort cannot contain is the thing tomato has.** Its duplication carriers read half
non-reference and its heterozygotes read half non-reference, and nothing else does. Real alignments
have a third population — positions where a fraction of the reads disagree in many accessions at once
because two parts of the genome map onto one place — and the trio showed the fit calling those
heterozygous with a posterior of 1.000 at a read share of a quarter
([`trio_heterozygosity_excess_2026-08-14.md`](trio_heterozygosity_excess_2026-08-14.md) §3). §5's
table is what that looks like from the cohort's side: the class's claim rate is highest exactly where
many accessions carry the allele.

---

## 7. Handing the class the depth reading, which is the obvious repair and is not one yet

The design has always named a second discriminator beside the cohort's genotype composition: **the
accession's own read depth around the position**, which is the only thing that works below about
twenty-five samples and the only thing at all at one sample. It has never been used on this cohort,
because the stored coverage summary keeps each accession's depth-against-GC curve but not each
window's own GC fraction, so the reading could not be corrected and would have been mostly GC
content. The audit computes the GC fraction from the reference during the walk, which is what
§10's change 6 would make unnecessary, and hands the fit `ln P(this depth | two copies) −
ln P(this depth | one)` at every position for every accession.

| | median accession's heterozygosity | less heterozygous than random mating predicts, median | the population's expected heterozygosity | the class's share of positions | passes |
|---|---:|---:|---:|---:|---:|
| the class off | 0.867 per kilobase | 0.789 | 4.895 per kilobase | — | 30 |
| the class on | **0.064** | 0.983 | 4.153 | 0.00789 | 30 |
| the class on, with depth | **0.404** | 0.888 | 4.559 | **0.03209** | 33 |

**Better, and wrong in a new way.** The depth reading returns the median accession from 0.064 to
0.404 heterozygous positions per kilobase — 47% of what the class-off fit books, where the cohort
pattern alone left 7%. But the class's share of positions goes the other way, from 1 position in 127
to **1 position in 31**: about 64,000 of the 1,999,404 kept positions, against a direct count of 24
to 468 per two million per accession.

**The two moves are one mechanism.** The reading enters the fit as a multiplier on the *carrier*
branch, so every position sitting in a window that reads high gets pushed toward the class whether or
not any accession's reads disagree with the reference there. The class swells to cover high-coverage
stretches in general; that raises the share it is fitted at, and — because the share is spread over
far more positions — it competes less hard for each individual heterozygous one. **The heterozygosity
comes back because the class is diluted, not because it has become selective.**

So the depth reading in this form is not the repair. What it would need is for the reading to bear on
*whether this position is in the class at all* rather than only on which accessions carry it, so that
a high-coverage window with no disagreeing reads anywhere is not a candidate.

---

## 8. What the class's default should be

**Off.**

The class was switched on by a measurement on drawn cohorts, where at 50 samples with no duplications
planted it invented a weight of 0.00003 and moved heterozygosity by half a per cent
([`contamination_floor_and_duplicated_class_2026-08-13.md`](contamination_floor_and_duplicated_class_2026-08-13.md)
§7). On 63 real accessions it books 0.00789 — 260 times that — and takes 93% of the heterozygosity
with it. §6 says the difference is not the shape of the frequency spectrum; §3 and §4 say the
positions it takes are not, in the main, duplications.

**What turning it off costs, and it is a real cost, not zero.** The class exists because the
two-class model has nowhere to put a position where a quarter to a half of the reads disagree in
everybody. On the benchmark human trio that is 59 positions in 449,489 carrying 79% of a
heterozygosity 26% above a real truth set
([`trio_heterozygosity_excess_2026-08-14.md`](trio_heterozygosity_excess_2026-08-14.md) §3). Off,
those positions go back to being counted as heterozygous. **The two errors are not the same size.**
Leaving them in puts the trio's heterozygosity **26% above** its truth set, and the drawn cohort of
§6.1 puts the same artefact at **11%** at tomato's own duplication rate; taking them out with this
class removes **93%** of tomato's. The class costs eight times what it saves.

**What would earn the default back.** A discriminator that separates a duplication from a
heterozygote at one accession at one position, measured on real reads rather than on a drawn cohort.
Depth is the candidate the design already names; §7 is what it does today, which is to recover half
the heterozygosity by claiming four times as many positions.

**And one thing that should not wait for that.** The class's per-position posterior now exists, so a
run that fits the class can say which positions it moved. **A heterozygosity that falls thirteenfold
should not be reportable without that list**, whatever the default becomes.

---

## 9. What this cannot say

- **There is no truth set for duplications in this cohort.** Every comparison here is between two
  ways of picking positions, not against truth. When §3 says the class books nine times as many
  positions as look duplicated, that is *nine times what a depth-and-read-share rule finds*, and the
  depth rule has its own misses: it is blind to a duplication shorter than the window it reads over,
  and to one whose two copies happen to agree at that base.
- **The depth reading has a floor of 0.7% to 3.4%**, which is the share of ordinary positions that
  read between 1.6 and 2.4 times normal. A claimed share of 14.7% is therefore about 12 points of
  signal over noise, and a claimed share of 9% on an accession whose floor is 3.4% is nearer 6.
- **Ten of the 63 accessions carry fewer than five reads a position**, and there the direct count of
  duplication-looking positions is restricted to the few positions with eight reads or more, so the
  25.5× ratio in §3's last row is an upper bound. The 8.9× on the six deepest accessions is the
  number to quote.
- **The rough SNP caller's 1.049 heterozygous positions per kilobase is not a truth set either.** It
  models neither duplications nor mismapping, so it is an upper bound with its own inflation. It is
  used here only to say that 0.064 is not in its neighbourhood.
- **The audit's near-half cell reads a position's own alternative-read share against its own depth**,
  taken as the middle of the range its stored code stands for. At the tomato panel's depths that
  range is one to three reads wide, so a position's share carries a little of the code's width.
- **The drawn cohort of §6.1 changed the frequency spectrum and nothing else.** Its duplication
  carriers still collect exactly twice the reads and still read exactly half non-reference, and every
  position it draws is one of three things. Real alignments are not so tidy — a duplication may be
  three copies, its two copies may differ at a rate the harness does not draw, and mismapping need not
  look like a raised error rate. So §6.1 rules out the spectrum's shape as the cause; it does not
  identify what tomato has that the harness lacks.
- **One drawn cohort per setting.** Differences of a percentage point there are not separable from the
  draw; the 11.3-point gap the class closes is.
- **One walk per arm of the real data.** The arms differ only in what the fit was told, so the
  differences between them are not a draw; the absolute values carry whatever the walk and the fit
  carry.

---

## 10. What I would change in the specifications

*Nothing under `spec/` or `arch/` was edited. These are the changes I would make, in the order I
would make them.*

**In [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md):**

1. **§2.2 — the class ships off, and the paragraph that turns it on needs the real-reads measurement
   beside the drawn one.** The section decides for the class on evidence that is entirely drawn or
   entirely about how many duplications exist; it has never carried what the class *does* on a real
   cohort large enough to identify it. On 63 tomato accessions it books 941 to 3,515 carrier
   positions per accession per two million where the depth reading finds 24 to 468, and takes the
   median accession's heterozygosity from 0.867 to 0.064 per kilobase.
2. **§2.2 — record that the cohort discriminator on its own does not select on depth, with the
   number that shows it.** Of the positions reading 35–65% non-reference at eight reads or more, the
   class takes 83% of those in a window at twice normal depth and 86% of those at ordinary depth. The
   document argues that the cohort's genotype composition "carries most of the benefit"; measured on
   real reads it carries no depth discrimination at all, and the two statements need to sit together.
3. **§2.2 — the per-position class posterior is an output of the fit, not an internal.** It is the
   only thing that turns the class's share into a list of positions, and every check in this report is
   a check on the list. Four bytes a position, beside the mismapped posterior §3.4 already keeps.
4. **§3.2 — the per-sample carrier rate belongs beside the heterozygosity in what the route
   produces.** `SampleGenotypeRates::duplicated_carrier` exists and nothing reports it; a run whose
   heterozygosity has fallen thirteenfold should say, in the same table, how many positions it moved
   and where they went.

5. **§2.2 — the depth discriminator's form is wrong, and the section should say what it currently
   does.** It enters the likelihood as a multiplier on the carrier branch, so a high-coverage window
   pushes a position toward the class whether or not any accession's reads disagree there. On the 63
   accessions it takes the class from 1 position in 127 to 1 in 31 while returning the median
   accession from 0.064 to 0.404 heterozygous positions per kilobase — half the heterozygosity back,
   bought by claiming four times the positions. The document names local relative coverage as the
   discriminator; it does not say where the reading is allowed to act, and that is the part that is
   wrong.

**In [`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md):**

6. **§4 — the coverage summary must keep each window's own GC fraction.** This is the change
   [`contamination_floor_and_duplicated_class_2026-08-13.md`](contamination_floor_and_duplicated_class_2026-08-13.md)
   §9 already proposed, and this report is the second use for it: without the window's GC the stored
   depth-against-GC curve cannot be looked up, so a relative depth reading is copy number times GC
   content. The audit here computes the GC fraction from the reference during the walk because the
   summary cannot supply it. One byte a window, against the one the mean depth already costs.
