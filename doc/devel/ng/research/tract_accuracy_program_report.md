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

*Sections P1–P2 and L1–L7 are opened as the program reaches them.*
