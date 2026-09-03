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

## P2 — where the never-offered truth sequences were lost

Status: **measured** — probe `benchmarks/ssr_hg002/src/never_offered_attribution.py`, 2026-09-03

**Pre-registration (written 2026-09-03, before the derivation is run).**

The superseded split of this class (434 under the old comparison: 268 no read carried it /
61 top-ploidy cut / 59 merge refused / 46 support bar) was produced by
`ng_tract_candidate_recall.py`, which carries two recorded defects
(`tract_genotype_accuracy_2026-09-03.md` §3.5): its window rule drops a truth insertion
anchored one base before the tract — where left-alignment puts every repeat-length gain —
and its truth reconstruction misses records left-aligned before the tract. This re-derivation
avoids both by construction: the missing truth sequences come from the corrected instrument's
own machinery (records collected within ten bases, reduced to minimal edits, left-aligned,
laid on the tract's reference), and the candidate side comes from the caller's own selection
code via `ng_candidate_selection_probe`'s per-tract dump (`NG_TRACT_DUMP`), which runs
`select_ssr` itself — every allele the merge tabled, whether some sample cleared the support
bar for it, and whether selection kept it.

targets: sizes L4 (the realigner's prize) and gates L7 (discovery's prize).
ceiling: not applicable — a measurement.

**The classes, each of the 463 never-offered tracts to exactly one:**

1. **merge refused** — no locus was built (`RefusedByMerge` in the dump);
2. **tabled and kept** — the truth sequence is among selection's survivors, yet the record
   cannot express it (§3.4c measured these as real disagreements, not artefacts; counted
   apart so they cannot inflate any lever's ceiling);
3. **top-ploidy cut** — in the table, cleared the support bar, dropped by selection.
   **This is the only class L7 (discovery) can reach**; its old count was 61 of 434;
4. **support bar** — in the table, never cleared the bar in any sample;
5. **never tabled** — no read carried the spelling into the merge. Subdivided by whether the
   truth's **length** was tabled (a spelling loss — the realigner's) or not (absent from the
   evidence or lost upstream). This class plus P1's unseen-in-raw is L4's pool.

**Controls:** the candidate dump is regenerated from the current binary and must reproduce
the Sep 2 stored dump (`tmp/attrib/tier_30x_candidates.tsv`) — selection code has not changed
since, so any difference is a red flag; and the join must find a dump entry for every one of
the 463 tracts, asserted, with the count of never-offered tracts derived by the attribution
equal to the verdict dump's.

**One extension, stated before it ran:** for the tracts whose truth length is absent from
the merge's table, `--bam-low` counts how many raw-aligned 30× spanning reads carry that
length, by P1's read rule — reads that carry it while the table does not are an admission or
realignment loss, the recall-side twin of P1's unseen-in-raw.

**Results (written after the runs they quote).**

Controls: the dump regenerated from the freshly built probe is **byte-identical** to the
Sep 2 stored one; the attribution derives exactly the verdict dump's 463 tracts (asserted);
every tract found its dump entry. Per-tract table `tmp/tract_program/p2_tracts.tsv`.

| class | homopolymer | period 2+ | total | old split (of 434) |
|---|---:|---:|---:|---:|
| merge refused | 7 | 7 | **14** | 59 |
| tabled and kept (real disagreements, §3.4c) | 5 | 1 | **6** | — |
| top-ploidy cut — **L7's whole reachable class** | 20 | 22 | **42** | 61 |
| support bar | 60 | 28 | **88** | 46 |
| never tabled, truth's **length** tabled | 42 | 46 | **88** | — |
| never tabled, length absent, **2+ raw reads carry it** | 65 | 49 | **114** | 67 ("alignment loss") |
| never tabled, length absent, 1 raw read | 7 | 7 | **14** | — |
| never tabled, length absent, 0 raw reads | 39 | 58 | **97** | 121 ("unrecoverable") |

**What this sizes:**

- **L4 (realigner), the program's largest prize by far.** 88 tracts where the table holds the
  right length in the wrong letters (the §6.2 corruption mechanism — an interruption or a
  leading base discarded), plus 114 where two or more raw reads carry a length the merge
  never tabled at all. That is **202 of the 463**, before adding P1's 74 spurious-side cases
  (65 unseen-in-raw + 9 spelling-only). The same caveat cuts both ways: "a raw read carries
  it" is the CIGAR-implied length, and the realigner legitimately re-spells some reads — the
  202 is the size of the disagreement surface between aligner and realigner, and L4's case
  reading decides how much of it is defect.
- **L7 (discovery): the gate number is 42**, down from the superseded 61 — every allele that
  cleared the support bar and was cut by the per-sample top-ploidy rule. Against it stands
  the unchanged risk side: the 225-strong spurious class it can only enlarge.
- **The support bar holds 88** — nearly double the old 46; the `hetfhom` sweep already
  measured that lowering it buys 2 tracts per 137 false candidates, so this stays closed.
- **97 tracts are a limit of the evidence** — no raw read at 30× carries the truth's length
  at all (the old 121, re-derived smaller).

---

## L4 — the realigner: what the reads are taken to say

Status: **diagnosed** — moved first in the order (owner's approval 2026-09-03); the fix
direction is a design decision now with the owner

**Pre-registration (written before any code is read or changed).**

This is a defect hunt in an existing component, not a build. The tract locus generator
realigns every read against the locus (`src/ng/locus_generation/ssr.rs`, using the unit-slip
whole-read aligner `src/ng/alignment/ssr_best_path_unit_slip.rs` — the bake-off's winner),
and the merge's allele table holds the spellings it produces. Three hand-verified cases show
those spellings corrupted (`tract_genotype_accuracy_2026-09-03.md` §6.2): a leading base
dropped (`chr3:33,877,690`), an interruption inside the tract discarded (`chr3:37,126,860`),
a one-read sequencing-error spelling kept over the true spelling carried by 12 reads of 14
(`chr11:37,147,255`).

targets: the aligner-vs-realigner disagreement surface Stage 0 measured —
  **202 never-offered tracts** (88 where the table holds the truth's length in the wrong
  letters + 114 where 2+ raw reads carry a length the table never held) and
  **74 spurious-side tracts** (65 where ng called a length no raw read spells + 9
  spelling-only differences), out of the 834 errors at 30×.
ceiling: **276 tracts**, and it will not be reached — the surface includes legitimate
  re-spellings (the realigner is allowed to disagree with BWA) and evidence limits; the
  ceiling is the size of the territory, not the expected win. No number smaller than the
  three verified cases is possible.
bar:
  - **default** (the fix ships) if a defect is found and fixing it flips more tracts right
    than wrong on both period classes at both depths, reported as verdict flips per rule 3,
    with the tomato behavioural gate run before it lands (a realigner change moves every
    run's records, not a parameter);
  - **discard** (the surface is legitimate disagreement) if case reading attributes the
    202+74 to correct re-spelling and evidence limits, with the count of read cases beside
    the verdict.

**Plan, as the program states it:** reproduce `chr3:33,877,690` in isolation — the reads
through the tract locus generator, watching where the sequence is lost — fix, then re-run
P2's attribution and the baseline pair. Any fix goes through the full plan-driven loop
(implement → review → apply fixes → commit).

**Diagnosis (written after the runs it quotes; the fix awaits the owner's design ruling).**

All three §6.2 cases are reproduced and understood, and they split into two mechanisms. Both
are junction events, and both defeat the delimiter the same way: the winning alignment
re-spells the read toward a pure motif run, destroying the non-motif evidence.

**Mechanism 1 — a real flank indel beside the junction (`chr3:33,877,690`).** The truth (and
BWA, and ng's own SNP path, which calls it correctly as `GT→G` het) put a one-base deletion
in the left flank's `TTT` and a `C` at the tract's first base. Reproduced end-to-end with a
one-interval run of the candidate probe (table holds `A×11`, `A×10`; 16 reads on the corrupt
spelling), then in a unit test against the exact 15-base window (kept `#[ignore]`d red as
`a_flank_indel_beside_the_junction_does_not_eat_the_tract_edge`). Two independent causes:

- the flank-side **junction guard** (7 columns here) makes the honest path's gap-open
  **unreachable** — a real indel within the guard window is inexpressible, the same failure
  the 4n flank ban was rejected for, recreated inside the guard's window;
- even unguarded, the cost structure prefers the corruption: a whole-unit slip open is
  ≈ −2.9 nats against the flank gap-open's ≈ −10.4, so "flank mismatch + tract contraction"
  beats "flank deletion + tract substitution" by ≈ 7.5 nats. The isolated-aligner test with
  no flank indel (`a_substitution_at_a_homopolymer_edge_stays_in_the_tract`, committed green)
  shows the substitution alone survives; it is the adjacent real indel that flips the path.

The end result double-counts the deletion (once as the SNP path's flank indel, once as a
tract contraction) and vanishes the `C` — it reaches no path at all.

**Mechanism 2 — the flank is itself a repeat (`chr3:37,126,860`, `chr11:37,147,255`).** Both
tracts are poly-A runs whose right neighbour is a sub-floor repeat the typing left as flank
(`(AAAG)×5`, `(GA)×5.5`). The truth's variant is one extra unit of that neighbour,
left-aligned into the tract span — inside the tract by the project's own convention. The
delimiter's "anchor" is five copies of a motif, so the inserted unit is absorbed into flank
matches for almost nothing, and the honest in-tract insertion is *also* inside the
tract-side guard window (these variants are junction events by construction). At
`chr11:37,147,255` the absorption even shifts the measured run: ng calls 13 A's on a 14-A
reference tract, homozygous.

**Sizes, on the baseline's 463 never-offered tracts** (probe run recorded in
`tmp/tract_program/l4_shapes.tsv`):

- **232 of 463** miss a truth sequence that is *not* a pure run of the tract's motif — the
  shape the delimiter cannot keep. Of the 88 right-length/wrong-letters tracts, 87 are this
  shape.
- The pool where the reads demonstrably carry what the table lacks: **117 tracts**
  (87 right-length + 30 length-absent with 2+ raw carriers, both interrupted), plus up to
  84 pure-length tracts with 2+ raw carriers where the delimiter re-measures the carried
  length, plus P1's 74 spurious-side cases.
- 110 of 463 have three or more reads carrying an input-alignment indel inside a 15-base
  flank — mechanism 1's reach.
- A crude "flank looks repeaty" flag does **not** separate errors from correct calls
  (0.60 among never-offered against 0.53 among right tracts): the discriminator is the
  variant engaging the junction, not the flank's texture.

**The decided arm (owner delegated the direction, 2026-09-03: "use your own criterium").**
Of the three directions — (a) spell the observation from the read's input alignment over the
tract span, (b) guard/pricing revision for corroborated flank indels, (c) compound-locus
typing — the build is **(a)**, chosen because it is the only one that reaches both mechanisms
at once, all three verified cases spell correctly under it, and it is judged by measurement
(verdict flips) rather than argument. (b) fixes mechanism 1 only; (c) is the structural
answer to mechanism 2 but moves the catalog, the ground and every benchmark, and stays a
follow-up informed by (a)'s residue.

**What (a) changes, precisely:** for a **complete** observation (the delimiter still rules
completeness, anchoring, the widen-and-retry, and the quality gate), the observation's bases
come from the read's own input alignment mapped over the tract's reference span, with the
settled junction conventions — an insertion at a boundary belongs to what follows, so a
left-junction insertion is the tract's and a right-junction insertion is the flank's. Where
the input alignment cannot serve (its aligned footprint does not bracket the tract plus one
base each side), the delimiter's own spelling stands, and partial observations are untouched.
Ordinary stutter reads spell identically under both (BWA left-aligns an in-tract indel to
the tract, so the span's content is unchanged); the two differ exactly where the delimiter
re-spells junction variation — the 232-tract shape.

**Bar for the arm:** verdict flips against the P0 baseline, net-positive on both period
classes at both depths, with the spurious-heterozygote class not enlarged beyond what the
fixed tracts justify; then the tomato behavioural gate before any default claim.

**Results (written after the runs they quote; arm `l4_input_spelling`, callsets
`tmp/tract_program/l4_arm/`).**

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| baseline | 0.8796 | 0.8692 | 0.8938 | 0.8780 |
| **input spelling** | **0.8981** | **0.8964** | **0.9114** | **0.9052** |
| gain (points) | +1.85 | +2.72 | +1.76 | +2.72 |

Errors fall from 834 to 679 at 30× and from 761 to 610 at 50×. **The flips, which are the
bar** (30×: 384 tracts changed verdict; join of the two `--verdicts-out` dumps):

- **197 tracts turn right**: 82 never-offered, 48 spurious heterozygotes, 48 wrong-some-other-way,
  13 collapsed heterozygotes, 6 previously silent or incomparable.
- **108 tracts leave right**: 59 become spurious heterozygotes (55 of them homopolymers),
  24 never-offered, 22 fall silent (no record at all), 3 other.
- Net by period class: homopolymer **+50**, period 2+ **+39**; at 50× **+51** and **+43**
  (sums match the arm's genotype-right deltas, +89 and +94). Both classes positive at both
  depths — the bar's first clause holds.
- The spurious-het class ends at 238 against the baseline's 225 (+13 net at 30×, +18 at
  50×) while total errors fall by 155 — within what the fixed tracts justify, and the 59
  created cases are recorded as L1's inheritance below.

**The counter-arm, built and rejected by measurement** (`l4_v2`,
`tmp/tract_program/l4_arm2/`): deferring reference-spelling inputs back to the delimiter — to
repair loci where the mapper places an ambiguous indel outside the tract
(`chr1:194,220,455`, where the same one-base contraction arrives spelled at three different
positions across reads) — nets only +50 tracts right against v1's +89, because at twice as
many loci the mapper's reference account is exactly right and the delimiter re-corrupts it:
at `chr1:104,262,346` every read carries a truth-confirmed flank deletion (`GC→G` hom), and
the deferral re-absorbs it as a spurious one-unit tract contraction — the original
double-count, resurrected. The deferral is reverted and the case pinned in the code comment.

**What the fix is, finally:** a complete observation's bases come from the read's input
alignment mapped over the tract span (junction conventions: an inserted run at a junction
belongs to the tract iff left-alignment cannot carry it across), taken only from reads whose
CIGAR carries an indel (an all-Match CIGAR is how a mapper spells a collapsed long allele —
the delimiter's one irreplaceable rescue, kept), and never in the widen-and-retry arm. The
delimiter still rules completeness, anchoring, quality, and partials. Validation:
6,034 lib tests, fmt, clippy all green in the container; the production-parity fixture and
the widen-recovery test forced the all-Match restriction and stand green.

**The residue, named so the next levers can claim it:**

- **59 tracts right→spurious-het** (30×): junk-shaped spellings near the true allele that
  the emission now misassigns — several are reads whose mapper account carries a stray
  junction mismatch bundled with the indel. This is L1's territory (junk mass near the
  called allele), and these tracts join its pre-registered target pool.
- **22 tracts right→silent**: junction insertions of foreign bases that the convention
  correctly expels from the tract but no path then owns (`chr1:163,385,956`: truth `T→TA`
  after the tract's last base; the baseline's tract record happened to carry it). The
  routing crack between the tract and generic paths at a region's first base — the mirror
  of main's `b6309954` — is a follow-up defect, not a spelling question.
- 26 never-offered tracts become silent rather than right — wrong either way, but silence
  is invisible to a record-level reader; same routing crack.

**Verdict: default-candidate.** The fix ships on this branch as behaviour (a defect fix,
not a parameter), and per rule 7 the *default* claim owes the tomato behavioural gate and
the owner's sign-off at Checkpoint 1, alongside the composed slate.

**P2 re-derived on the fixed caller** (probe rebuilt, dump regenerated, control asserted;
`tmp/tract_program/l4_p2_tracts.tsv`): the never-offered class is 369 (was 463). The
never-tabled pool falls 313 → 203; tracts where 2+ raw reads carry a length the table lacks
fall **114 → 30**; right-length-wrong-letters falls 88 → 63 — the fix consumed over half the
realigner surface, and 100 of the remaining length-absent tracts have no raw carrier at all
(the evidence limit, stable). **The top-ploidy cut rises 42 → 61**: truth alleles that now
survive into the table are being cut by the per-sample ploidy rung — discovery's (L7's)
reachable class grew back to its old size, and its gate is re-armed with that number.

---

*Sections L1–L3 and L5–L7 are opened as the program reaches them.*
