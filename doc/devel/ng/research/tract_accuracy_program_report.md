# The tract-accuracy program's report — one section per lever

**The plan is [`tract_accuracy_program.md`](tract_accuracy_program.md)**; this file accumulates
its results. A section is opened with its pre-registration *before* the lever is built, and its
measurements and verdict are written *after* the runs they quote — never during.

Every number is HG002, GIAB tandem-repeat benchmark, 30× unless stated, scored by
`benchmarks/lib/tract_qual_experiment.py` (corrected comparison, self-test green), genotype
accuracy letter-for-letter on the tract's own bases.

**Section template** (copy for each new lever):

```
## L<n> — <name>
Status: pre-registered | built | measured | VERDICT: default / optional / discard
Pre-registration (written before building):
  targets: <error class and count>
  ceiling: <the arithmetic maximum it could reach>
  bar:     <what earns which verdict>
Arms and results: <tables, with arm labels and run paths>
Verdict: <one of three, with the number that decides it, and what it still owes>
```

---

## Sections closed before the program started

These levers were measured while the comparison was being corrected (2026-09-03, sweep arms
under `tmp/slip_sweep/`, `tmp/conc_sweep/`; tables in
[`ng_tract_genotype_improvement_2026-09-02.md`](../../reports/ng_tract_genotype_improvement_2026-09-02.md) §2
and commit `a753c0b0`). Recorded here so the program does not re-open them.

### junk-strength — the outlier weight's value

**VERDICT: default, kept at 0.05.** Across 0.01–0.30 (thirty-fold) homopolymer accuracy moves
inside 0.4 points with the 0.05–0.30 stretch flat within 0.1; the classes trade (141 spurious
heterozygotes at 0.01 → 108 at 0.30, collapsed 37 → 56). 0.05 takes the whole period-2+ gain
and two thirds of the homopolymer one. **Still owed:** the joint re-sweep with L1's shape
change, and the tomato range gate — a floor under every read behaves differently at 3 reads
than at 30, and nothing has tested it there.

### flat slippage settings — slip share, shorter share, fall-off, hand-set length gradient

**VERDICT: discard, all four.** Slip share over 0.02–0.40: accuracy moves half a point while
the two error classes swing 2× and 9× in opposite directions; shipped 0.10 sits beside the peak
(0.05, ahead by 0.14 points at period 2+ and behind by 0.09 at homopolymers — inside the trade,
not above it). Shorter share and fall-off: every setting within 0.1 points of shipped. The best
of three hand-set length-rising shapes: +0.04.

### the prior's length spectrum — fitting it

**VERDICT: discard.** The class it would re-balance splits 110 homozygous-reference against
132 homozygous-non-reference (probe `tmp/levers/spurious_het_shape.py`, 30×): a spectrum peaked
at the reference length fixes the first bucket and worsens the second, a wash by construction.
The truth's spectrum is real (79 chromosomes in 100 at the reference length against a flat
shape's 11 in 100) — it is informative about the genome and not a lever for this error.

### the prior's concentration — the one het/hom dial the tract path has

**VERDICT: default, kept at 1.0.** The tract prior is a marginalised Dirichlet–multinomial;
with K candidates and total belief C, the heterozygous share of prior mass is
(K−1)·C / (K·(C+1)) — 42% at K=6, C=1. Swept 0.10–8.0 at the shipped outlier weight
(`tmp/conc_sweep/`): homopolymer peak is exactly the shipped 1.0; period 2+ gains 0.04 points
at C=8 for 19 more spurious heterozygotes. Third dial to trade the same two errors, third to be
found already at its optimum.

### base quality on the tract path

**VERDICT: discard (owner, 2026-09-03).** A read's own base qualities are never read at a
tract — the emission scores letters against a per-stratum fitted substitution rate. The signal
is real (mean base quality 35.3 inside tracts against 34.0 outside, `readmodel`
investigation) but carrying per-read qualities into the tract emission is an architectural
change the owner rules out.

### a stricter candidate bar, a GQ floor, an allele-balance collapse rule

**VERDICT: discard**, from the `hetfhom` investigation — each costs more correct calls than it
removes wrong ones (best case: 37 of 137 removed for 617 true alleles lost; GQ 30 withdraws
413 correct calls per 43 wrong ones and converts errors to no-calls, not to right calls).
Counted under the old comparison; re-opened only if a program lever changes the landscape they
were measured on.

---

## P0 — the baseline pair and the verdict dump

Status: **done** — all four bars met, 2026-09-03

**Pre-registration (written 2026-09-03, before the instrument edit).**

P0 is not a lever; it is the measuring stick every lever is read against. Two deliverables:

1. **`--verdicts-out` on the instrument** — one row per tract the truth calls, carrying the
   tract's coordinates, its period class, and which counter it landed in (right, no-call,
   not comparable, truth allele never offered, spurious heterozygote, collapsed heterozygote,
   wrong some other way), plus whether the repeat lengths alone were right. Rule 3 (verdict
   flips, never headlines alone) needs this on every arm; today the per-tract comparison lives
   only in throwaway scripts under `tmp/`.
2. **A fresh `--defaults` baseline at 30× and 50×** — the stored callsets under
   `benchmarks/ssr_hg002/results/ng/` were run when the outlier default was 0.01; the shipped
   default moved to 0.05 (commit `073b678c`) and the program's baseline must be the shipped
   caller, freshly run.

targets: no error class — the instrument itself.
ceiling: not applicable.
bar (all four must hold before anything downstream is scored):

- the self-test is green before and after the edit;
- the edit leaves all three existing tables **byte-identical** on a re-score of the stored
  30× callset (the harness-change control, rule 1);
- the fresh 30× `--defaults` run agrees with the sweep's existing `outlier0.05` arm
  (`tmp/slip_sweep/outlier0.05/calls.vcf`) — the same setting reached two ways: same genotype
  row, and **zero verdict flips** between the two callsets in the new dump;
- the 30× headline reproduces the corrected baseline, 0.8796 homopolymer / 0.8692 period 2+.

**Results (written after the runs they quote).**

The instrument now writes `--verdicts-out`: one row per truth-called tract —
coordinates, period class, verdict (`right` / `no_records` / `no_call` / `not_comparable` /
`never_offered` / `spurious_het` / `collapsed_het` / `wrong_other`), and whether the repeat
lengths alone matched. `benchmarks/lib/tract_verdict_flips.py` joins two dumps on the tract
and prints the flip crosstab. All four pre-registered bars:

1. **Self-test**: green before and after the edit (7 pins; two new checks assert the verdict
   rows themselves).
2. **Harness-change control**: the stored 30× callset re-scored before and after the edit —
   `calibration.tsv`, `sweep.tsv`, `genotype.tsv` all **byte-identical**
   (`tmp/tract_program/control_pre/` against `control_post/`).
3. **Same setting, two roads**: the fresh 30× `--defaults` run against the sweep's
   `outlier0.05` arm — the VCF **data lines are identical** and the flip join reports
   **6,993 tracts, 0 flipped**. This also confirms the freshly rebuilt binary reproduces the
   callset the sweep's arm was scored from.
4. **The headline reproduces**: 0.8796 / 0.8692 at 30×, to the fourth decimal.

**The program's fixed baseline** (fresh `--defaults`, shipped outlier weight 0.05; callsets
`tmp/tract_program/baseline/HG002_{30x,50x}.raw.vcf`, verdicts
`tmp/tract_program/verdicts.tsv`, arm label `baseline`):

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| sequence accuracy | 0.8796 | 0.8692 | 0.8938 | 0.8780 |
| repeat-length accuracy | 0.8829 | 0.8749 | 0.8969 | 0.8843 |
| never offered | 245 | 218 | 215 | 198 |
| spurious het (called het, truth hom) | 129 | 96 | 124 | 98 |
| collapsed het (called hom, truth het) | 40 | 20 | 23 | 13 |
| wrong some other way | 51 | 35 | 51 | 39 |

834 errors at 30×, 761 at 50×. The spurious-heterozygote class holds at the new default what
it showed at 0.01: **225 at 30× against 222 at 50×** — 1 case in 100 fewer where total errors
fall by 9 in 100 — so the class the levers target is still the one more reads do not buy down.

Flip check run on the way (0.01 stored callset → 0.05 arm, `tract_verdict_flips.py`):
28 tracts flipped, 17 spurious hets fixed against 4 collapsed hets created — the trade §1 of
the plan describes, now visible tract by tract.

---

## P1 — what the spurious-allele reads have in common

Status: **measured** — probe `benchmarks/ssr_hg002/src/spurious_read_provenance.py`, 2026-09-03

**Pre-registration (written 2026-09-03, before the probe is built).**

This is the measurement that aims the program: at each of the baseline's **225
spurious-heterozygote tracts at 30×** (verdict `spurious_het` in P0's dump), pull the reads
that carry the spurious length and ask whether they are independent evidence, and whether the
locus keeps producing that length when there are ten times the reads.

targets: no error class — it partitions the 225 across the levers.
ceiling: not applicable.

**Definitions, fixed before anything is measured:**

- *The spurious length.* The called pair's sequence(s) absent from the truth's pair, taken
  over the tract's own bases; its length in bases. A spurious sequence with the **same**
  length as the truth's (a spelling difference only) goes to its own bucket — reads binned by
  length cannot arbitrate it.
- *A spanning read.* Primary, mapped, covering the tract plus one base each side.
- *A read's tract length.* Aligned bases placed on the tract's own positions, plus inserted
  bases anchored from the base before the tract through the tract's second-to-last position —
  the same rule at both depths, so any edge convention cancels in the comparison.
- *Clustered* (the reads share an origin): among the k spurious-length reads, either
  (a) all k sit on one strand and the chance of that under the tract's own strand mix is
  below 1 in 20, or (b) at least half of them share one exact (start, strand, template
  length) signature with k ≥ 2 — the shape of a PCR duplicate family, which these BAMs do
  not flag (checked: `samtools flagstat` reports 0 duplicates at both depths).
- *Persists at 300×.* At least 3 reads of the 300× alignment carry the spurious length AND
  its share there is at least **half** its 30× share. The 30× BAM is a subsample of the 300×
  one, so the 30× reads alone reproduce only ~a tenth of their 30× share — the bar demands
  genuinely new reads.

**The partition, each tract to exactly one bucket:**

| bucket | criterion | points at |
|---|---|---|
| spelling-only | spurious length = truth length | the realigner (L4), not a read-count lever |
| clustered | clustered, not persistent | L2 — n reads that are not n pieces of evidence |
| locus-real | persistent, not clustered | L3/L4 — the locus really yields the length |
| both | clustered AND persistent | listed one by one (expected rare) |
| sampling noise | neither | no lever — an unlucky draw more reads dilute |

**Controls:** the probe's own tract set must equal P0's verdict dump — 225 tracts, asserted
in the probe, not eyeballed. Beside the partition, the probe re-reports the share-of-reads
distribution so it can be read against the 0.01-callset's measured shape (158 of 242 at
3 reads in 10 or more).

**What re-ranks the levers:** a clustered majority sends L2 up; a locus-real majority
confirms L3/L4; a sampling-noise majority says the class is priced about right and the
program's weight shifts to the 463 never-offered errors (P2, L4).

**Results (written after the runs they quote).**

**One deviation from the pre-registration, made after first sight of the data and applied
before any verdict:** 65 tracts have **zero** 30× reads carrying the spurious length at all,
which the persistence bar did not anticipate — a share of zero makes "at least half the 30×
share" vacuously true, and the read questions (strand, family, persistence *of the carriers*)
are unanswerable with no carriers. Those tracts get their own bucket, `unseen_in_raw`, now
encoded in the probe itself. Control passed: the probe derives exactly the dump's 225
`spurious_het` tracts, asserted in code.

**The partition of the 225** (per-tract table `tmp/tract_program/p1_tracts.tsv`, whole cases
in `tmp/tract_program/p1_cases.txt`):

| bucket | homopolymer | period 2+ | total | what it says |
|---|---:|---:|---:|---|
| clustered | 0 | 0 | **0** | **L2's premise is false** — not one tract's spurious reads share a strand beyond their tract's own mix or an identical (start, strand, template-length) signature |
| locus-real | 60 | 73 | **133** | the reads genuinely carry the length, and it holds at 300× |
| unseen-in-raw | 52 | 13 | **65** | no raw-aligned read spells the called length at 30×; for 36 of the 65, none of a median 210 reads at 300× does either |
| sampling noise | 11 | 7 | **18** | median share falls from 0.062 at 30× to 0.019 at 300× |
| spelling-only | 6 | 3 | **9** | same length, different letters — reads binned by length cannot arbitrate |

**The locus-real 133 are not stutter-shaped, and that is the finding that re-aims the
program.** Their spurious share is median **0.458 at 30× and 0.486 at 300×** (105 of 133 hold
0.30 or more at 300×); the carriers sit on both strands at independent start positions; their
MAPQ is a uniform 70 (checked at two cases; no paralog pile-up); 99 of 133 sit exactly one
motif unit from the truth, mostly at 10–24-unit tracts. Half the reads carrying a second
length at 300× is not slippage — the measured per-read slip rate tops out at 9 in 100 at
30-unit homopolymers — it is the reads and the assembly-based truth genuinely disagreeing
about the sample. Hand-checked whole: at `chr1:77568472` (21-base poly-A) the truth writes
`CAA→C` on **both** haplotypes (homozygous 19), and 110 of 241 reads at 300× spell 20 bases.
Whether that is an assembly error on one haplotype or a real property of the DNA is a
geneticist's question, not a measurement's; what the measurement rules out is any read-weighing
lever fixing these without also refusing genuine one-unit heterozygotes that look identical.

**The unseen-in-raw 65 point at the realigner (L4), and give it a second charge.** The plan's
L4 section anticipated this: "P1 may show it also manufactures spurious second lengths." A
called allele that no raw alignment spells even once in a median 210 reads at 300× (36 tracts) is
either the realigner spelling reads differently than the aligner does — legitimately or not —
or an allele manufactured in candidate construction. This is the same machinery as §6.2's
three hand-verified corruption cases, from the other side.

**The share distribution, for continuity with the 0.01-callset shape** (measured there from
ng's own AD: 158 of 242 at 3 reads in 10 or more): by raw-read share, 107 of the 216
length-distinct tracts are at 3 in 10 or more. The two measures differ — AD counts reads ng's
realigner assigns, this counts reads whose raw alignment spells the length — and the 65
unseen-in-raw tracts are the bulk of the gap.

**How this re-ranks the levers:**

- **L2 (read-independence): premise false, measured.** Zero clustered tracts of 225. Per its
  own pre-registration rule in the plan, L2 stops at arm (a), run cheaply for the record.
- **L3 (stutter-em): ceiling shrinks.** The one-unit signature it targets is real (99 of 133)
  but at a median share of 0.46 — far beyond what a pulled-back per-locus re-fit can or
  should explain. The honest remaining target is the 28 locus-real tracts whose 300× share
  is under 0.30.
- **L4 (realigner): stock rises.** 65 unseen-in-raw + 9 spelling-only = 74 of the 225 now
  sit in its territory, before P2 re-derives its share of the 463 never-offered.
- **L1 (junk-shape): target unchanged but small.** The thin-evidence tail it absorbs is the
  18 sampling-noise tracts plus part of the 86 wrong-some-other-way.
- **The 133 locus-real are, on this benchmark, a wall** — reads and truth disagree; no
  caller lever reaches them without collateral damage. They cap what any composition of
  L1–L3 can show on the spurious-het class.

*(the geneticist's read on the locus-real class is Checkpoint 0's question to the owner)*

---

*Section P2 and L1–L7 are opened as the program reaches them.*
