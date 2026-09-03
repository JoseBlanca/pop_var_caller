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

*Sections P0–P2 and L1–L7 are opened as the program reaches them.*
