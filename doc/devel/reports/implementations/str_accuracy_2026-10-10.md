# STR-segment accuracy against GATK: six fixes measured on the stuttering benchmark

*Implementation report, 2026-10-10, branch `str-accuracy`, merged to main. Benchmark:
[`str_stutter_benchmark.md`](../../ng/spec/str_stutter_benchmark.md), first pass
[`str_stutter_benchmark_first_pass_2026-10-10.md`](str_stutter_benchmark_first_pass_2026-10-10.md).
Every caller run here used the simulator's true stutter written into the parameters, one group per
library (no fit: the fit currently pools every library into one group, a defect being fixed
elsewhere). Scratch, not committed: `tmp/acc/` in the branch's worktree.*

---

## 1. What was measured

The same 16 simulated runs as the first pass: PCR and PCR-free libraries, 1 and 50 samples, 5, 15,
30 and 60 reads a position. The unit is one sample at one STR segment (a *pair*), scored right,
wrong, missed or false call by `scripts/ng_str_stutter_score.py`. The 262 hexanucleotide segments
of motif AAAAAC, AAAAAG or AAAAAT are left out throughout, as in the first pass, so each run has
2,738 segments: 1,125 variant and 1,613 reference pairs at one sample, 56,621 and 80,279 at 50.

Each change was built and run on its own, against the build before it, and committed only after
its measurement. The one change to how reads are measured was also run on GIAB HG002's
tandem-repeat benchmark (36,497 truth records; 30× and 50×), scored by
`benchmarks/lib/tract_qual_experiment.py`.

## 2. Where main lost, at one sample

Of the pairs GATK called right and main did not (PCR, one sample, 15×: 110 against 25 the other
way), the causes were, in order of size: a called heterozygote dropped because its site quality was
under 1 (44); a deletion running from the repeat into its flank (22); the sample's allele never
made a candidate (23); a wrong genotype written (21).

## 3. The changes

| commit | change | spec | one sample, 8 runs: right / false calls added | 50 samples, 4 runs (5×, 15×): right / false added |
|---|---|---|---|---|
| 3919b1a7 | a tract's site quality uses the tract's own prior, not the SNP/indel one | `calling_quality.md` §8.1 | +167 / +34 | +2 / 0 |
| d82ad1a9 | one read off the motif grid can no longer refuse a tract | `candidate_alleles_ssr.md` §7 | +7 / 0 | 0 / 0 |
| a2dfa4f2 | a repeat length the reads earned keeps its best spelling when errors split them | `candidate_alleles_ssr.md` §5 | +63 / 0 | +151 / +7 |
| 2cf69b8d | a candidate with no whole repeat is dropped, not the tract | `candidate_alleles_ssr.md` §5 | 0 / 0 | +122 / 0 |
| c9800222 | a read whose alignment never reaches a tract is no observation of it | `locus_generation_ssr.md` §3 | 0 / 0 | +83 / −6 |
| c0b55125 | a read carrying a deletion across a tract's edge crossed the tract | `locus_generation_ssr.md` §3, `candidate_alleles_ssr.md` §7 | +76 / 0 | +1,122 / +16 |

The last also gained 138 right genotypes on HG002 over 30× and 50×, with one tract a depth going
from right to wrong. HG002's allele-exact comparison counts about 75 more false positives a depth
under it: 70 of the 114 new records at 30× sit at tracts whose genotype became right, written as
two records (flank piece, tract piece) where the truth writes one.

**Not adopted:** letting one crossing read make a candidate instead of two. At 5× it added about as
many false calls as right genotypes (PCR, one sample: 216 newly right, 252 newly false), and
genotype quality separated them only partly.

## 4. Against GATK, before and after

Right variant pairs / false calls. *Before* is main at 4b528afe as scored in the first pass, with
the true stutter.

| run | before | now | GATK |
|---|---|---|---|
| PCR, 1 sample, 5× | 437 / 2 | 506 / 14 | 529 / 9 |
| PCR, 1 sample, 15× | 940 / 0 | 1,022 / 17 | 1,025 / 18 |
| PCR, 1 sample, 30× | 1,080 / 0 | 1,106 / 3 | 1,095 / 21 |
| PCR, 1 sample, 60× | 1,101 / 0 | 1,114 / 1 | 1,107 / 23 |
| PCR-free, 1 sample, 5× | 490 / 2 | 544 / 2 | 541 / 2 |
| PCR-free, 1 sample, 15× | 984 / 0 | 1,029 / 1 | 1,044 / 0 |
| PCR-free, 1 sample, 30× | 1,092 / 0 | 1,106 / 0 | 1,112 / 1 |
| PCR-free, 1 sample, 60× | 1,102 / 0 | 1,112 / 0 | 1,120 / 0 |
| PCR, 50 samples, 5× | 38,926 / 1,432 | 39,745 / 1,449 | 25,546 / 2,004 |
| PCR, 50 samples, 15× | 52,787 / 336 | 53,216 / 338 | 49,873 / 7,309 |
| PCR, 50 samples, 30× | 55,265 / 76 | 55,739 / 76 | 53,267 / 9,445 |
| PCR, 50 samples, 60× | 55,945 / 14 | 56,354 / 14 | 53,927 / 9,961 |
| PCR-free, 50 samples, 5× | 40,888 / 844 | 41,522 / 853 | 27,522 / 172 |
| PCR-free, 50 samples, 15× | 54,010 / 81 | 54,347 / 81 | 52,603 / 546 |
| PCR-free, 50 samples, 30× | 55,846 / 11 | 56,198 / 11 | 55,784 / 404 |
| PCR-free, 50 samples, 60× | 56,035 / 0 | 56,433 / 0 | 56,281 / 229 |

Per run, one sample: 1,125 variant and 1,613 reference pairs; 50 samples: 56,621 and 80,279.

## 5. What is left

- **One sample, low depth.** Of the pairs GATK still gets and main does not, nearly all have 3 or
  fewer reads crossing the segment for one of the sample's alleles; the evidence is in reads that
  end inside the segment, which score genotypes but cannot propose an allele, and which can only
  prove a longer allele, not a shorter one.
- **Deletions across a segment's edge.** At one sample main gets 16–17 of the 27 such pairs right
  at 15× and above, GATK 26–27. At 50 samples (PCR, 15×) 347 of the 1,085 carrier pairs are still
  not right, 289 of them at the 18 of 62 carried edge-deletion alleles whose sequence the records
  still cannot spell (before the rule: 33 of 62, 563 pairs); `scripts/ng_str_accuracy/edge_deletions.py`
  counts it. One proposal, not built: let a segment's read
  spelling take in a few flank bases each side, so such a deletion is one allele of one record.
- **Representation.** An edge deletion is now written as two records; allele-exact comparisons
  count that as a miss and a false positive.
