# The stuttering benchmark, first pass: ng against main, GATK and HipSTR on STR segments of known genotype

*Implementation report, 2026-10-10, branch `tract-ceiling`. Spec:
[`str_stutter_benchmark.md`](../../ng/spec/str_stutter_benchmark.md). What was built: the simulator
(`examples/ng_str_stutter_benchmark.rs`), the driver (`scripts/ng_str_stutter_benchmark.sh`), the
scorer (`scripts/ng_str_stutter_score.py`), the grid (`scripts/ng_str_stutter_first_pass.sh`) and
its summary (`scripts/ng_str_stutter_grid_summary.py`). Outputs under
`tmp/stutter_bench/first_pass/` (not committed): each run's `score_cells.tsv`, and the grid's
`grid_cells.tsv`, `grid_vcfeval.tsv` and `headline.tsv`.*

---

*Outcome, owner, 2026-10-10: on segments with no variant beside or across their edge — where the
boundary rule of §2.2 never applies — the branch is behind main in all 8 runs at 50 samples (PCR 15×,
fitted: 92.43% against 93.64%) and ahead only at one sample and 5–15×. So calling both ways is not
kept: branch `tract-ceiling` is not merged, and its spec `tract_dual_path.md` stays there. The
compound-motif fix, the segmenter speed fixes and this benchmark's tools were moved to main.*

---

## 1. What was measured

A simulated cohort whose STR genotypes we drew, read through libraries that stutter at rates taken
from the literature and this project's own fits (spec §3), mapped with minimap2 and called
by every caller. **Sixteen runs**: PCR-free and PCR libraries, 1 and 50 samples, 5, 15, 30 and 60
reads a position. Each run: 4 Mb of genome, 3,000 STR segments (500 per motif length, from the
calling floor to 30 or 15 copies), SNPs and indels between them.

**The unit scored is one sample at one STR segment** (a *cell*): the caller's records over the
segment and 30 bases either side are applied to the reference, and the cell is **right** when they
spell the sample's two true sequences — the spec's rule, edits slid left and counted per copy, so
an STR record and generic records saying the same thing score alike. Other outcomes: **wrong**
(a variant called with another sequence), **missed** (called homozygous reference, or left
uncalled), and **false call** (a variant called where the sample carries none). At 50 samples a
run has 150,000 cells, of which 62,127 carry a variant; at one sample, 3,000 and 1,241.

**Callers**, each called once per run:

| name here | what it is |
|---|---|
| branch, defaults | this branch, `call-from-alignments --defaults` |
| branch, fitted | this branch, `call-from-psps` with the fitted parameters |
| branch, true stutter | this branch, the defaults with the simulator's own stutter written in (no fit) |
| main, defaults / fitted / true stutter | main (4b528afe) the same three ways, with the branch's same two parameters files |
| GATK | HaplotypeCaller per sample (GVCF), CombineGVCFs, GenotypeGVCFs |
| HipSTR | v0.6.2, all samples jointly over the segments; at its defaults it skips a segment with fewer than 100 reads over the cohort |
| HipSTR, fixed stutter | the same with `--def-stutter-model --min-reads 1`, its own answer for small cohorts |

Main against the branch is the comparison that measures **calling each STR segment both ways**:
main calls each segment by the STR path alone, fits with the same code, and reads the branch's
parameters files unchanged (its echoed file equals the fit line for line). One difference is not
that change: main lacks today's compound-motif fix (db7388ec) and types the 262 hexanucleotide
segments of motif AAAAAC, AAAAAG or AAAAAT as ordinary sequence. **Every number below leaves those
262 out**, so both builds are scored on the same 2,738 segments (56,621 variant cells at 50
samples, 1,125 at one).

**The fit was computed once per cohort**, at 30×, and reused at 5×, 15× and 60× (owner): a cohort is
the same samples at every depth. So at 5× the fitted runs use parameters fitted on 30× reads, which
a 5× cohort could not have fitted itself.

## 2. What it shows

### 2.1 At one sample, calling both ways gets more right than main

Variant cells called right, one sample (1,125 cells):

| | PCR 5× | 15× | 30× | 60× | PCR-free 5× | 15× | 30× | 60× |
|---|---|---|---|---|---|---|---|---|
| branch, true stutter | **44.7%** | **88.2%** | **97.9%** | **98.8%** | **47.4%** | **90.8%** | **98.7%** | **99.7%** |
| main, true stutter | 38.8% | 83.6% | 96.0% | 97.9% | 43.6% | 87.5% | 97.1% | 98.0% |
| branch, fitted | 40.1% | 82.0% | 95.7% | 98.4% | 43.4% | 86.9% | 97.4% | 99.4% |
| main, fitted | 33.7% | 79.5% | 95.0% | 97.7% | 40.2% | 85.8% | 97.0% | 98.0% |

The branch is ahead of main at every depth with every parameter set, by 1 to 6 points. **It pays in
false calls**: 16 against 2 at PCR 5× with the true stutter, 15 against 0 at 15×, out of 1,875
reference cells. At 60× main makes about 22 wrong genotypes where the branch makes 2 to 8.

### 2.2 At 50 samples, the branch withholds calls beside segments and comes out behind main

Variant cells called right, 50 samples (56,621 cells):

| | PCR 5× | 15× | 30× | 60× | PCR-free 5× | 15× | 30× | 60× |
|---|---|---|---|---|---|---|---|---|
| main, true stutter | 68.75% | **93.23%** | 97.61% | **98.81%** | 72.21% | 95.39% | **98.63%** | 98.97% |
| branch, true stutter | **68.91%** | 92.62% | 96.74% | 97.66% | **72.57%** | 95.29% | 98.62% | **99.42%** |
| branch, true stutter, placeholders as absent | 69.16% | 93.15% | 97.63% | 98.64% | 72.63% | 95.48% | 98.87% | 99.60% |
| main, fitted | 68.31% | 92.74% | 97.21% | 98.47% | 72.01% | 95.21% | 98.62% | 98.97% |
| branch, fitted | 67.68% | 91.67% | 96.28% | 97.55% | 71.31% | 94.48% | 98.45% | 99.36% |
| branch, fitted, placeholders as absent | 67.90% | 92.15% | 97.06% | 98.43% | 71.33% | 94.63% | 98.62% | 99.54% |

**The cause is one rule.** Where an SNP/indel edit crosses a segment's edge and the stutter weight is
0.95 or more, the branch writes the STR answer and marks the flank bases that edit touched as
missing for every sample — records with no alternative allele, `FILTER=strBoundary`, every sample
`./.` (`tract_dual_path.md` §8). The scorer counts a cell under such a record as uncalled, which is
what the VCF says. In PCR runs these placeholders are few — 30 to 80 records a run — but each covers
all 50 samples:

| PCR, 50 samples, fitted | 5× | 15× | 30× | 60× |
|---|---|---|---|---|
| placeholder records | 30 | 39 | 54 | 60 |
| branch: variant cells uncalled | 1,042 | 673 | 877 | 961 |
| main: variant cells uncalled | 1,020 | 127 | 135 | 13 |

Some of what they withhold is real: at the cells checked by hand, main had called a SNP in that
flank which the sample carries (chr1:721709 A→T, chr3:317079 C→T), and the branch replaced it with
missing data. Others withhold nothing the sample carries — its STR genotype is right, and three
flank bases where it carries nothing are written as missing.

**With the placeholders treated as absent records, the branch and main are level** (true stutter,
PCR: 69.16% against 68.75% at 5×, 93.15 against 93.23 at 15×, 97.63 against 97.61 at 30×, 98.64
against 98.81 at 60×). The two also err differently: the branch makes fewer wrong genotypes (PCR 60×,
true stutter: 481 against 630) and misses more variant cells (289 against 46).

In PCR-free runs the placeholders are rarer (4 to 16 records a run) and the branch is level with or
ahead of main as written.

**At default parameters the branch is ahead of main at both cohort sizes** (50 samples, PCR: 67.0%
against 65.8% at 5×, 96.4 against 96.1 at 30×) — the both-ways combination helps most where the
stutter model is wrong.

### 2.3 The true stutter beats the fit everywhere

The run with the simulator's true stutter and no fit calls more right than the fitted run in every
one of the 16 runs, by 0.1 to 6.1 points (largest at one sample: PCR 1 sample 15×, 88.2% against
82.0%), and makes fewer false calls in 7 of the 8 runs at 50 samples (PCR 15×: 418 against 952;
the exception, PCR-free 5×, 1,007 against 968). **So on
simulated data the fit can be skipped** (owner's criterion): writing the true stutter in costs
seconds where the fit costs 2.7 to 3.2 hours at 50 samples, and is the better answer. It also says
the fit leaves accuracy on the table — but the fit code here is the one being replaced in a
parallel effort, so that is a note for that work, not for this branch.

### 2.4 Against GATK and HipSTR

Variant cells right / false calls among the 87,873 reference cells, 50 samples (all 3,000 segments
here would change these by under a point; the 2,738 shown):

| | PCR 5× | 15× | 30× | 60× | PCR-free 5× | 30× |
|---|---|---|---|---|---|---|
| branch, true stutter | 68.9 / 1,516 | 92.6 / 418 | 96.7 / 89 | 97.7 / 16 | 72.6 / 1,007 | 98.6 / 29 |
| GATK | 45.1 / 2,004 | 88.1 / 7,309 | 94.1 / 9,445 | 95.2 / 9,961 | 48.6 / 172 | 98.5 / 404 |
| HipSTR | 63.0 / 5,191 | 89.0 / 2,602 | 92.8 / 1,483 | 93.6 / 1,134 | 67.5 / 1,894 | 93.9 / 108 |
| HipSTR, fixed stutter | 62.7 / 7,307 | 87.0 / 5,079 | 91.8 / 4,050 | 92.8 / 3,449 | 66.3 / 2,005 | 93.5 / 100 |

- **ng calls the most cells right at 50 samples at every depth and in both designs.**
- **GATK's false calls grow with depth in PCR libraries** (2,004 at 5×, 9,961 at 60×), and sit on
  long mono- and dinucleotides at the stutter plateau (at 30×, 8,027 of its 9,445 at 16 repeats or
  more). Checked by hand: it calls stutter reads as an allele — a homozygous-reference sample at a
  run of 16 T's called heterozygous for a one-base deletion from 9 of 31 reads, at GQ 99.
- **At one sample GATK is level with or ahead of ng at 15× and above** (PCR 15×: 91.1% against
  88.2%; PCR-free 30×: 98.8 against 98.7), with more false calls in PCR (18 against 15 at 15×).
- **HipSTR at its defaults calls nothing at one sample** (fewer than 100 reads at every segment).
  With its fixed stutter model it is ahead of ng at one sample and 5× (PCR 51.6% against 44.7%) at
  the cost of 316 false calls against 16.

### 2.5 Away from the segments nothing changed

vcfeval outside the segments' windows: the branch and main give identical true positives, false
positives and false negatives in every one of the 16 runs (e.g. PCR 50 samples 15×: 57,535 / 101 /
164 both). The change is confined to STR segments, as designed.

## 3. Cost

| run | total | of which the fit | GATK | HipSTR |
|---|---|---|---|---|
| PCR, 50 samples, 30× | 4 h 33 min | 3 h 09 min | 19 min | 10 min |
| PCR-free, 50 samples, 30× | 3 h 51 min | 2 h 42 min | 13 min | 11 min |
| PCR, 50 samples, 60× (fit reused) | 2 h 10 min | — | 24 min | 20 min |
| PCR, 1 sample, 30× | 16 min | 14 min | 46 s | 3 s |

The grid took 14 hours, about 6 of them fitting; main's pass, 2.5 hours more. Without the fit — §2.3
— a 50-sample run at 30× is about 1.5 hours.

## 4. What this does not show

- **One seed, one genome, one cohort per design and size.** Differences of a few tenths of a point
  at 50 samples (a few hundred cells) are within what another seed could move; the one-sample
  differences of 1 to 6 points are 10 to 70 cells out of 1,125.
- **Clusters of repeats, satellites, interruptions and real sequence context are absent** (spec §7).
  The simulated genome has no paralogs, so the duplication filter's effect is only its false drops.
- **The heavy-PCR, clean and one-outlier designs and 10-sample cohorts are not run yet** (spec §4.1's
  full grid). The one-outlier design is the case spec §1 flags.
- **The placeholder count is a scoring choice as well as a calling one.** The scorer counts a
  placeholder as uncalled; §2.2 shows both readings.
- **ng's `--explain-loci` output was kept** in every run but not yet read against the truth (spec
  §6's check of whether the stutter weight favoured the right model at each disagreement).

## 5. For the owner

1. **The `strBoundary` rule costs the branch its lead at 50 samples.** Without it, the branch is level
   with main at 50 samples in PCR and ahead in PCR-free; with it, 0.5 to 1.2 points behind in PCR at
   15× and above. The rule is the spec's (§8: a boundary-crossing edit at weight ≥ 0.95 is written
   from the STR answer, the rest missing), and this was one of the three numbers the benchmark was
   to test. **Recommendation: write the generic call instead of missing data over the flank bases**
   — keep the STR record, and keep the SNP/indel records the generic path called beside it — and
   re-run the PCR 50-sample cohort to confirm. Decision is the owner's: it is a spec change.
2. **Skip the fit on simulated runs**: use the true stutter (§2.3).
3. **The branch's extra false calls at one sample** (16 against main's 2 at PCR 5×, out of 1,875
   reference cells) are small in count; worth a look at the next pass, not now.
