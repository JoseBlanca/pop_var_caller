# `--call-satellites` against GATK on GIAB HG002, over satellite ground

**2026-10-18, branch `call-satellites`.** The accuracy check of plan
[`call_satellites.md`](../../ng/impl_plan/call_satellites.md), step 6: what `--call-satellites`
calls on satellites, and how much of it agrees with the GIAB truth, beside GATK.

**One sample, HG002, at three depths.** Nothing here measures the cohort machinery; at one sample
the option's work is all in the SNP/indel model, which is what is being judged.

## The answer

**On satellite ground `--call-satellites` finds about as many true variants as GATK, a little
less precisely, and gets their genotypes right less often.** Over the Tandem-Repeat benchmark's
satellite ground at 30 reads a position, counting a variant found when its alleles match the truth
whatever the genotype:

| | SNPs found (of 1,940) | SNP calls that are true | other variants found (of 513) | other calls that are true |
|---|---:|---:|---:|---:|
| GATK | 793 (41%) | 95% | 278 (54%) | 65% |
| pvc `--call-satellites` | 740 (38%) | 91% | 229 (45%) | 57% |
| pvc default | 0 | — | 1 | 1 of 3 |

At 10 reads a position the two find the same share of SNPs (24% each) and GATK stays ahead on the
rest (39% against 32%). Requiring the genotype too, pvc gets 383 SNPs right at 30× against GATK's
542, and 291 against 372 at 10×.

**Both callers miss most of what the truth holds on satellites**, at every depth: GATK finds 41%
of the SNPs at 30×. Short reads do not place well in long arrays, and the misses are shared.

**The limit the spec warned about barely binds.** Loci refused for being wider than
`--max-locus-span` (100 bp): 3 over the whole 30× run with the option, against 2 without; none at
10×. pvc finds 45% of the non-SNP variants on satellites, against GATK's 54% — so the satellites'
length variation is not mostly lost to width, as `command_line.md` §4.4.1 expected; it is missed
about as often as GATK misses it, for the reason both miss SNPs.

## The data

| set | reads | ground | truth | satellite ground in it |
|---|---|---|---|---|
| 300× | HG002, 1,000 regions, 5.08 Mb (`benchmarks/human_genome_bottle/`) | GIAB small-variant benchmark regions | GIAB small-variant truth over those regions | 12 satellites, 1,464 bases; 13 truth variants |
| Tandem-Repeat 30× and 10× | HG002 subsampled from 288× (`benchmarks/ssr_hg002/`) | GIAB Tandem-Repeat v1.0.1, 50,000 Tier-1 regions, 6.09 Mb | that benchmark's assembly-based truth | 328 satellites, 51,459 bases; 2,439 truth variants (1,959 SNPs) |

**Satellite ground** is what `pvc dev type-regions` types `satellite` under the default settings,
clipped to the benchmark regions. The runs agree with it: each `--call-satellites` run's report
says it handed over exactly those bases (1,464; 51,459).

**The 300× set is too small to judge by**: 13 truth variants. There, GATK finds all 13 with no
false call, pvc `--call-satellites` all 13 with one false call, and the default run 1. The small
benchmark regions avoid long repeats, which is why the Tandem-Repeat set carries the result.

**GATK** is 4.6 HaplotypeCaller with `-stand-call-conf 0`: the existing 300× VCF, and for the
Tandem-Repeat set a run over the satellite intervals padded by 150 bp at each depth.

## How it was scored, and why two ways

Two methods, beside each other, on the same files:

- **`vcfeval`** (RTG Tools 3.13) replays both call sets onto the reference and compares sequences,
  in genotype mode and allele mode (`--squash-ploidy`). The tables above are its numbers.
- **The project's exact matcher** (`benchmarks/lib/compare_to_truth.sh`'s method: restrict to the
  BED, left-align, split multi-allelic records, then require position and alleles to agree).

**On satellites the exact matcher undercounts both callers by about three times**, and that is the
case for not using it there: at 30× it finds 236 of GATK's true SNPs where `vcfeval` finds 793, and
219 of pvc's where `vcfeval` finds 740. Both callers write neighbouring substitutions in an array as
one multi-base record, which the matcher's normalisation does not split, so a correct call counts as
a false call and several misses.

| 30× SNPs, true calls | exact matcher | `vcfeval`, alleles | `vcfeval`, genotypes |
|---|---:|---:|---:|
| GATK | 236 | 793 | 542 |
| pvc `--call-satellites` | 219 | 740 | 383 |

## Every number

`scores.tsv`, written by the scripts below; *tp* is truth variants matched (for `calls`, the count
of call records on satellite ground), *fp* call records unmatched, *fn* truth variants unmatched.

| set | caller | measure | class | tp | fp | fn |
|---|---|---|---|---:|---:|---:|
| 30× | GATK | calls | SNP / other | 288 / 669 | | |
| 30× | GATK | vcfeval, genotypes | SNP | 542 | 72 | 1,313 |
| 30× | GATK | vcfeval, genotypes | other | 252 | 240 | 251 |
| 30× | GATK | vcfeval, alleles | SNP | 793 | 38 | 1,147 |
| 30× | GATK | vcfeval, alleles | other | 278 | 151 | 235 |
| 30× | pvc `--call-satellites` | calls | SNP / other | 256 / 435 | | |
| 30× | pvc `--call-satellites` | vcfeval, genotypes | SNP | 383 | 125 | 1,557 |
| 30× | pvc `--call-satellites` | vcfeval, genotypes | other | 191 | 300 | 322 |
| 30× | pvc `--call-satellites` | vcfeval, alleles | SNP | 740 | 77 | 1,200 |
| 30× | pvc `--call-satellites` | vcfeval, alleles | other | 229 | 170 | 284 |
| 30× | pvc default | calls | SNP / other | 0 / 3 | | |
| 10× | GATK | calls | SNP / other | 176 / 461 | | |
| 10× | GATK | vcfeval, genotypes | SNP | 372 | 38 | 1,568 |
| 10× | GATK | vcfeval, genotypes | other | 162 | 188 | 351 |
| 10× | GATK | vcfeval, alleles | SNP | 472 | 20 | 1,468 |
| 10× | GATK | vcfeval, alleles | other | 198 | 113 | 315 |
| 10× | pvc `--call-satellites` | calls | SNP / other | 223 / 334 | | |
| 10× | pvc `--call-satellites` | vcfeval, genotypes | SNP | 291 | 102 | 1,636 |
| 10× | pvc `--call-satellites` | vcfeval, genotypes | other | 115 | 234 | 397 |
| 10× | pvc `--call-satellites` | vcfeval, alleles | SNP | 479 | 67 | 1,461 |
| 10× | pvc `--call-satellites` | vcfeval, alleles | other | 163 | 139 | 350 |
| 10× | pvc default | calls | SNP / other | 0 / 2 | | |

`vcfeval` counts 1,940 truth SNPs where the exact matcher counts 1,959: it sets aside truth
records it cannot evaluate inside the regions.

**Where pvc's genotypes differ from the truth, nothing stands out.** Of the 180 calls at 30×, of
every class, whose alleles match and genotype does not, 69 are called heterozygous, 40 homozygous and 71
with two alternative alleles; the truth has 1,941 heterozygous and 478 homozygous variants there.

## Reproducing it

Scripts in `scripts/call_satellites_accuracy/`, paths written for this worktree:

1. `pvc dev type-regions --reference <GRCh38 with its catalogue beside it> --regions <BED>` for
   each set's regions, then `satellite_bed.py <typed.tsv> <regions.bed> <out.bed>`, then the BED
   sorted into the reference's contig order.
2. `run_pvc.sh` (host): pvc with and without `--call-satellites`, `--defaults`, 6 threads.
3. `run_gatk.sh` (dev container, benchmarks mounted with `DEV_EXTRA_MOUNT`).
4. `score.sh` (dev container): both methods, into `scores.tsv`. RTG's `format` built the
   reference's SDF first.

Run on the branch after merging `main` (48bd2360); the scores were identical before the merge.
