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

## L1 — junk-shape: where the junk mass sits

Status: pre-registered

**Pre-registration (written 2026-09-03, after L4 landed and before any L1 code).**

**The baseline for this lever is the L4-fixed caller** (arm `l4_input_spelling`):
0.8981 / 0.8964 at 30×, 679 errors — 369 never-offered, 238 spurious het, 44 collapsed het,
28 wrong-other.

**The geneticist's criterion, from the owner (2026-09-03):** much of the junk is expected to
be **hidden duplications** — single-copy in the reference, multi-copy in the sample — where a
mutation in one sample-copy manufactures an artefactual heterozygote, at a distance that can
be **more than one motif unit**. Two measured facts frame it:

- The refreshed P1 partition on the L4 baseline (probe re-run, control asserted, 238 tracts):
  **184 locus-real** (persistent at 300×, both strands, MAPQ 70), 22 unseen-in-raw,
  16 sampling-noise, 12 spelling-only, 4 clustered-and-persistent. The locus-real class shows
  **no coverage excess** at 300× (median 175 spanning reads against 182 over all scored
  tracts) — the classic collapsed-duplication signature (roughly doubled depth) is absent, so
  if these are duplications they are ones whose extra copy's reads are split or lost, not
  piled on.
- **The stray tail is multi-step, as the owner predicted**: of the 52 spurious alleles
  carried by under 2 reads in 10, only 12 sit one unit from the truth — the spread reaches 36
  units, with most inside 7. A junk shape peaked at ±1 unit is contraindicated; the contest
  is uniform (shipped) against a moderate decay concentrating floor mass within a few units.

targets: the **52 stray-tail spurious heterozygotes** plus part of the 28 wrong-other. The
  184 well-supported locus-real cases are explicitly NOT claimed — a junk term absorbs stray
  reads, and a second length carried by half the reads at 300× is not stray.
ceiling: **80 tracts** (52 + 28), and it will not be reached — some strays are real slippage
  the stutter model owns.
bar: verdict flips against the L4 baseline, net-positive on both period classes at both
  depths; the collapsed-het class not enlarged beyond what fixed strays justify (a stronger
  junk floor is a het-suppressor by construction); λ re-swept jointly with any shape change
  (rule 6); the simulator sanity run (planted junk absorbed without moving clean loci); and
  the tomato behavioural gate before any default claim — a floor under every read binds
  hardest at 3 reads a position, which HG002 cannot see.

**Arms, in order:** (b) the junk weight split per period class, seeded from the measured
junk rates (1 read in 2,300 at homopolymers, 1 in 209 at period 2+ — the opposite ordering
to what the single 0.05 implies) and swept around them; (a) the uniform `U` replaced by a
geometric decay with distance in motif units from the nearest candidate, decay swept from
steep to nearly flat with uniform as the control; (c) the winner of each combined, λ
re-swept on top.

*(measurements to follow the runs — nothing below this line was written before them)*

**Arm (b), measured — and closed without code.** Each locus is scored under its own period
class independently, so per-class optima read directly off a global-λ sweep; the split needs
plumbing only if the peaks separate. They do not: re-swept on the L4 baseline
(`tmp/tract_program/l1_sweep/`, six values 0.01–0.30 at 30×, the two contenders at 50×), the
homopolymer optimum is **0.20** and period 2+ is flat from 0.05 to 0.20 (within 0.04 points),
so a global 0.20 captures everything a per-class pair would. **L4 moved this dial's optimum**:
pre-L4 the 0.05–0.30 stretch was flat within 0.1 points; post-L4, 0.05 → 0.20 is worth
**+0.32 / +0.04** at 30× (flips: 25 tracts right — 24 of them spurious hets — against 11 new
collapsed hets and 6 lost records; net positive both classes) and **+0.58 / +0.16** at 50×
with **no collapse cost there** (17 collapsed homopolymer hets, identical to baseline).
The λ=0.05 control row reproduces the L4 baseline to the fourth decimal. Recommendation to
the slate: **default λ = 0.20**, pending the tomato behavioural gate — a floor under every
read binds hardest at 3 reads a position, which HG002 cannot see.

**Arm (a), measured — the shape is a weaker copy of the strength dial.** The knob
(`repeat_tract_junk_decay_per_unit`, a stated constant defaulting to 1.0 = the shipped
uniform) is built end-to-end and proven inert at its default twice over: the row takes the
pre-change expression verbatim at 1.0 (bit-identity asserted in a unit test), and a fresh
`--defaults` run under the decay build is **byte-identical** to the pre-decay L4 arm
(`CONTROL-BYTE-IDENTICAL`, `tmp/tract_program/l1_decay_run.sh`). The sweep
(decay 0.85/0.70/0.50/0.30 × λ 0.05/0.20 at 30×, `tmp/tract_program/l1_decay/`):

- at λ 0.05, decay 0.50 is the shape's best: 0.9008 / 0.8964, **flips +6 net** (22 right
  against 16 leaving) — strictly inside what uniform λ 0.20 buys (+8 net, 0.9013 / 0.8968)
  on the same spurious-for-collapsed trade;
- on top of λ 0.20 the shape only destroys: its mildest setting nets **−6** (32 right
  against 38 leaving), and the steepest reaches 0.8917 with 88 collapsed homopolymer hets.

**Verdict for the shape: discard as a default** — dominated at every grid cell by the
uniform with the right λ, which is what the owner's multi-step-junk reading predicts: junk
lands at many distances, so a uniform floor was already approximately the right shape. The
built knob ships inert at 1.0; whether it is removed outright or kept as the fallback point
on the trade curve is decided in the Checkpoint 1 slate — its one conceivable use is if
λ = 0.20 fails the 3× tomato gate, since a decayed 0.05 reaches similar HG002 numbers with
less floor mass per read. The simulator sanity run is waived for a discarded shape and owed
before any future adoption.

**Verdict for the weight (the owner's junk-strength lever, §5's owed joint re-sweep):
default-candidate λ = 0.20**, replacing 0.05 — +0.32 / +0.04 at 30×, +0.58 / +0.16 at 50×,
no 50× collapse cost — pending the tomato behavioural gate and the owner's slate.

**Set aside with the owner (2026-09-03): a per-locus λ re-fit in the EM loop.** Proposed,
examined, and ruled out before building: its reachable band is the 5–20%-share strays only
(the pull-back that makes it safe also makes it unable to touch the 46%-share wall), its
failure mode is un-calling real heterozygotes, and at 3 reads a position it is inert by
design. The 184-tract wall's real discriminator is the **cohort** — heterozygote excess
across samples, repeated off-ratio allele balance, Hardy–Weinberg violation — signals a
single-sample benchmark cannot exercise; that is a future cohort-level lever outside this
program, connected to the existing hidden-paralog census work.

---

## L2 — read-independence: n agreeing reads are not n pieces of evidence

Status: **premise measured false — discard proposed to the owner (their lever), for the
Checkpoint 1 sitting; nothing built**

The lever family (the identical-observation discount, the freebayes-style aggregate factor,
the GATK-style per-read cap, the beta-binomial) exists to stop reads that share an origin
from counting as independent evidence. P1 measured the premise directly, twice: on the P0
baseline, **0 of 225** spurious-heterozygote tracts show a strand or duplicate-family
clustering signal; re-measured on the L4 baseline, **0 of 238** show clustering alone and
only 4 show clustering beside persistence. The reads behind the spurious class sit on both
strands at independent start positions — they are exactly as independent as real evidence,
so a discount can only reweigh good and bad evidence alike, which is the strength dial λ
already does with one number. The caps' own pre-registered ceilings concur: the aggregate
factor is bounded by a tenth of the unexplained mass, and a per-read floor bounds one read
where these classes are many.

**The number for the ruling: 0 clustered of 238, with 4 ambiguous.** If the owner rules the
discard, the family is closed unbuilt; if not, arm (a) (the one-line discount, swept) is the
cheapest measurement to run first.

---

## L3 — stutter-em: the locus's own slippage

Status: pre-registered

**Pre-registration (written 2026-09-03, before the re-fit body is built).**

The calling loop's per-locus slippage re-fit is designed and refused-at-runtime today
(`calling_em_loop.md` §5.1; `SlippageRefitConfig`, 50 pseudo-counts, 20 slipped reads,
`SlippageRefitNotBuilt`). This lever builds the body inside the existing interface and sweeps
the three designed pull-back settings. **Scope ruling (owner, 2026-09-03): slippage only** —
the per-locus λ re-fit that could share this machinery was examined and set aside (L1's
section records why).

targets: the locus-real spurious heterozygotes whose 300× share is **under 0.30** — the band
  a pulled-back re-fit can plausibly explain as locus-specific slippage. On the L4 baseline:
  **60 tracts** (46 homopolymer, 14 period 2+), 31 of them exactly one unit off. The
  ≥0.30-share wall (128 tracts) is explicitly not claimed — the pull-back cannot and should
  not reach it.
ceiling: **60**, and part of the 44 collapsed heterozygotes if a locus fitted *below* its
  stratum re-arms real one-unit hets.
bar: flips against the L4-plus-λ0.20 baseline once the slate fixes λ (or against L4-λ0.05
  with λ re-swept jointly if the slate is still open when L3 runs — rule 6 either way);
  net-positive both period classes both depths; the simulator run (its slippage is settable,
  so the re-fit must recover a planted per-locus rate on it — the mechanism check rule 5
  requires); the 242-tract subset scored directly per the plan; and the interaction with
  fitted per-stratum rows (L6) measured at stage end.

**Results (written after the runs they quote; arms under `tmp/tract_program/l3_arms/`).**

The re-fit body is built inside the designed interface (`slippage_refit.rs`, the round driver
in `summarise_condition.rs`; production's formulas with the pull-backs as settings, the
spec's posterior-weighted attribution in place of production's hard assignment; granularity
mirrored from `em.rs` — one pooled count set per locus, one level multiplier, per-cell shape
pull-back — with the deviations documented at the module head). 6,060 lib tests, fmt, clippy
green. **Two identity controls hold**: frozen (`rounds = 0`) is byte-identical to the L4 arm,
and rounds switched on at `--defaults` is byte-identical too — every cell is a shipped
constant there, outside the re-fit, so the rounds measure and adopt nothing. The re-fit
therefore engages only above fitted per-stratum rows, which folds the planned L3×L6
interaction into the measurement by construction.

On the fitted rows (`fitted_slippage_hg002_30x.toml`), 30×:

| arm (λ 0.05) | homopolymer | period 2+ | spurious het | collapsed het |
|---|---:|---:|---:|---:|
| fitted rows, frozen | 0.8932 | 0.8983 | 175 + 94 | 21 + 7 |
| + re-fit, designed pull-backs (50/20, ≤3 rounds) | 0.8934 | 0.8986 | 172 + 93 | 23 + 7 |
| + re-fit, free (HipSTR's zero pull-back) | 0.8950 | 0.8927 | 135 + 88 | 46 + 26 |

**The pulled-back re-fit is a near-no-op: 7 verdict flips of 6,993, net +1 tract**
(4 spurious hets fixed against 2 new collapsed and 1 silenced); at λ 0.20 the same
(+0.02 / 0.00). The free setting is the EM degeneracy the pull-back exists to prevent,
measured: collapsed hets double and period 2+ loses 0.56 points.

**Verdict proposed: discard as a default** — the owner's lever, so the number goes to the
Checkpoint 1 sitting: **+1 net tract of 6,993 at the designed settings**. The machinery stays
built and frozen at zero rounds (building it was the spec's own requirement — §12's Q2 is now
answered), with the env switch remaining experiment-only.

**A finding that re-aims L6, recorded here because this measurement produced it:** the fitted
per-stratum rows themselves — worth +0.18 before L4 — now **cost homopolymers 0.4–0.5
points** against plain defaults (64 flips: 28 right→spurious at homopolymers against
16 collapsed fixed and 10 silent tracts recovered; net ≈ −4). Two readings, unresolved: the
corrected observations no longer carry the corruption the low fitted rates were fitted on —
**the rows are stale, measured on the pre-L4 caller's observations** — or low per-stratum
rates genuinely over-trust one-unit reads. L6's real question is now "re-fit the strata on
the fixed caller and re-measure", not "build the fit-mode command around the old rows".

---

## The tomato behavioural gate (rule 7), for the λ = 0.20 candidate

**Run 2026-09-03** on the 63-accession tomato cohort (the `tomato1` bench slice, ~2 Mb, ~3
reads a position), both arms on the current binary, one cohort invocation each, the 0.20 arm
replaying the 0.05 run's own parameters file with only the outlier weight edited
(`tmp/tract_program/tomato_gate/`). No truth set — what moved:

| | λ = 0.05 | λ = 0.20 |
|---|---:|---:|
| generic records | 227,881 | 227,881 — **identical**, as the dial must leave them |
| tract records | 966 | 938 (−28, 3 in 100) |
| het share among tract sample-genotypes | 0.068 | 0.062 |
| hom-alt / hom-ref / no-call | 0.094 / 0.752 / 0.086 | 0.099 / 0.752 / 0.088 |
| tract QUAL quartiles | 247 / 568 / 1,585 | 198 / 453 / 1,181 |

**No breakage shape**: no no-call surge, no record collapse, no het wipe-out at 3 reads. The
cost side for the owner's standing conservatism question: about **1 tract heterozygote in 11
is no longer called het** at 3× (0.068 → 0.062), and tract QUALs sit about a fifth lower.

---

## Checkpoint 1 — the owner's rulings (2026-09-03)

1. **λ = 0.20 adopted as the shipped default** (gate passed; the 1-in-11 tract-het trade at
   3× accepted).
2. **The junk-decay knob removed** — it gains nothing over the weight, so the parameter goes.
3. **L2 (read-independence / the two caps): explanation requested before ruling** — see the
   chat record; the measured premise stands at 0 clustered of 238.
4. **L3 (stutter-em): disabled, code kept** — rounds ship at zero, machinery stays.
5. **L5/L6: re-measure first** — re-fit the per-stratum rates on the fixed caller, decide
   then.
6. **L7 (new-alleles): approved to proceed.**

Main was merged back in under the levers; its parallel junction-ownership work includes the
routing-crack fix L4's residue named (`328a1a2b` — the flank claims the insertion the tract
refuses). Measured on the post-adoption baseline: **13 flips against the λ-0.20 arm, every
one a silent tract turned right** (6 homopolymer, 7 period 2+), nothing else moved — 13 of
the 22 tracts L4's residue named, recovered by the parallel fix with zero breakage.

**The program's baseline, restated after Checkpoint 1's adoptions** (fresh `--defaults`,
λ = 0.20, decay knob removed, main's junction fix in; arm `adopted`,
`tmp/tract_program/adopted/`):

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| program start (P0) | 0.8796 | 0.8692 | 0.8938 | 0.8780 |
| **after Checkpoint 1** | **0.9015** | **0.8970** | **0.9173** | **0.9071** |
| gain (points) | +2.19 | +2.78 | +2.35 | +2.91 |

Errors at 30×: 834 → **665** (369 never-offered, 214 spurious het, 56 collapsed het,
26 other). The tomato gate stands as run (its λ arm is now the default).

---

## L7 — new-alleles: discovery, wired and measured

Status: pre-registered

**Pre-registration (written 2026-09-03, before the wiring is built).**

The decision half is built and tested (`calling/inference/discovery.rs`, 12 tests, E1); the
wiring belongs inside `select_ssr` per E1's findings: discovery is a **pre-pass, not a round
wrapped around the loop** — the tract locus generator already realigns every read, so the
eligible set is a function of the observations and the candidate table alone, and a second
round over the same evidence admits nothing (asserted by an E1 test). The spec's §4.1 premise
("look against the converged posteriors") is contradicted by that and is amended rather than
the code bent to it — recorded at E1 and standing.

targets: the **61 top-ploidy-cut tracts** (P2 re-derived on the fixed caller) — a truth
  allele cleared the support bar and the per-sample ploidy rung cut it; discovery admits it
  back when 2+ reads and 15%+ of the sample's spanning reads carry it.
ceiling: **61**, minus whatever the admission bar refuses.
risk, named: every admitted length can only enlarge the spurious-het class (214 under the
  adopted λ 0.20); the bar is that flips stay net-positive on both period classes at both
  depths with the spurious class's growth smaller than the never-offered class's shrinkage.

**Results (written after the runs they quote; arm `l7_discovery`, `tmp/tract_program/l7_arm/`).**

The wiring landed as E1's findings direct: the pre-pass runs inside `select_ssr` over the
merge's own support rows (a support row *is* a complete realigned observation), per sample
against the shipped bar (2 reads AND 15% of that sample's spanning reads, never pooled across
samples), admitted alleles face the shared cap and truncation rules, and a new
`DiscoveryMode::BeforeTheLoop` names the setting while the two posterior-round modes stay
refused. Off is proven byte-identical (`CONTROL-OFF-IDENTICAL`); 6,061 lib tests, fmt, clippy
green; switch `NG_TRACT_DISCOVERY=1` (parameters-file plumbing owed on a keep).

Measured against the adopted baseline:

| | 30× | 50× |
|---|---|---|
| accuracy (hom / p2+) | 0.9026 / 0.8974 (+0.11 / +0.04) | 0.9186 / 0.9074 (+0.13 / +0.03) |
| flips | 31: **11 right-gains** (10 never-offered, 1 silent) against 4 leaving right | 21: **10 right-gains** against 2 leaving (both to spurious) |
| never-offered | −23 (369 → 346) | −17 |
| spurious het | **unchanged** (214) | +2 |
| wrong-other | +19 — errors converted from never-offered into offered-but-mis-genotyped | +7 |

**The bar is met on every clause**: net-positive both period classes both depths (hom +5 /
p2+ +2 at 30×; +6 / +2 at 50×), and the spurious class's growth (0 and 2) is far under the
never-offered shrinkage (23 and 17). The pre-registered ceiling was 61; discovery reached 25
of those tracts outright (11 + 14 converted) — the rest stay cut by the bar or mis-genotyped
over the enlarged set, which is the next lever family's territory, not this one's.

**Verdict: default-candidate, pending the owner's sign-off — the tomato gate passed.**
On the 63-accession cohort at ~3 reads (`tmp/tract_program/tomato_gate_l7/`), discovery on
against off changes **2 records of 228,852** and nothing else, byte for byte: at
`SL4.0ch06:22,887,980` and `SL4.0ch08:30,899,713` one interruption-carrying spelling is
admitted with coherent genotype and AD shifts (a `1/1` becomes `2/2` on the fuller spelling
with the extra allele's read in AD). At three reads the 2-reads-and-15% bar almost never
fires beyond what the merge already tabled — the mechanism is inert exactly where thin
evidence makes admission dangerous, which is what the range commitment asks of a
depth-dependent lever.
