# Improving ng's genotype accuracy at repeat tracts — the research program

**Date:** 2026-09-03. **Asked for by the owner**: try every lever we have collected, one at a
time — implement it, measure it, and then decide whether it becomes the shipped default, an
option a run can ask for, or is discarded with the number that discards it. The assistant runs
the program as independently as it can, and consults the owner only at the points §7 names:
a shipped default changing, one of the owner's own levers being discarded, or a question that
needs a geneticist's judgement rather than a measurement.

**The companion report is
[`tract_accuracy_program_report.md`](tract_accuracy_program_report.md)** — one section per lever
attempted, written after its runs and never before them. This file is the plan and stays stable;
that file accumulates.

**Context a reader needs first:**
[`tract_genotype_accuracy_2026-09-03.md`](tract_genotype_accuracy_2026-09-03.md) — how the
numbers are produced, and §3 for why this quantity has been measured wrong five times.
Every figure below is from the corrected comparison (its §3.4b) unless it says otherwise.

---

## 1. What we are trying to move, and from where

**The baseline is ng's shipped defaults** — the stutter model at HipSTR's constants, the
repeat-tract outlier weight at 0.05, one chromosome's worth of flat prior belief over each
tract's candidate lengths, nothing fitted. On GIAB's HG002 tandem-repeat benchmark at 30×,
genotype accuracy scored letter-for-letter on the tract's own bases:

| | homopolymer | period 2+ |
|---|---:|---:|
| sequence accuracy | 0.8796 | 0.8692 |
| repeat-length accuracy | 0.8829 | 0.8749 |
| tracts scored | 3,861 | 2,822 |

**The 834 errors, by what would have to change to fix them:**

| | count | what reaches it |
|---|---:|---|
| a truth sequence was never offered as a candidate | 463 | the realigner and candidate selection — nothing in the calling loop |
| called heterozygous, truth homozygous | 225 | the calling loop — this is the program's main target |
| called homozygous, truth heterozygous | 60 | discovery would; most levers trade against it |
| wrong some other way over a set holding the right alleles | 86 | mixed |

**The wall the program starts at.** Three settings — the slip share, the outlier weight, and the
prior's concentration — each trade a spurious heterozygote against a collapsed one, and each
already sits at or beside its own optimum on this benchmark. Re-weighting the two answers
against each other is exhausted. What the spurious heterozygotes actually are (measured on the
outlier-0.01 callset, 242 cases; the shipped 0.05 leaves 225 of the same shape):

- in **158 of 242** the weaker of ng's two alleles carries **3 reads in 10 or more** — a
  well-supported second length, not a stray read;
- **138 of 242** sit **exactly one motif unit** from the truth;
- the class **does not shrink with depth** (242 at 30×, 236 at 50×, while total errors fall
  from 852 to 774) — an error more reads do not buy down;
- at homopolymers the rate climbs with tract length: 28 in 1,000 at 10–14 motif units,
  73 in 1,000 at 25 or more.

So the surviving levers are the ones that change **what the reads are taken to say** — how much
evidence n agreeing reads constitute, what the stutter model expects at this locus, where the
junk mass sits, and what the realigner reports — not how two answers are weighed.

**And the range caveat stands over everything**: every number above is one human sample at
30× and 50×. The caller is committed to one-sample-to-thousands and 3× to hundreds
(`design_principles.md` §0), and §2's rules say what that costs each experiment.

---

## 2. The rules every experiment follows

1. **One instrument.** `benchmarks/lib/tract_qual_experiment.py`, self-test green before any
   scoring; runs driven by `sweep_tract_parameters.sh` (one arm) or
   `RESCORE_ONLY=1 run_tract_qual_experiment.sh` (re-derive everything). A harness change is
   followed by the byte-identical control run before any arm is trusted.
2. **Pre-register before building.** Each lever's report section is opened *before* its code is
   written, with three things: the error class it targets, the **arithmetic ceiling** (how many
   of the 834 it could possibly reach if it worked perfectly), and the success bar. A mechanism
   that cannot produce its number is caught here, not after a week of building.
3. **Verdict flips, not headlines.** Every arm reports fixes and breaks against baseline —
   which tracts changed verdict and which way — beside the headline. The program's own history
   says why: the one proposed correction judged by its headline alone raised it 4.7 points and
   corrected nothing. The instrument's per-tract verdict dump serves this (P0).
4. **Both period classes, both depths, always.** A lever that buys homopolymers by selling
   period 2+ is reported as exactly that, never netted.
5. **The simulator is a sanity check, not a verdict.** It has no mapping error, no
   interruptions and no aligner losses, and it overstates stutter levers by an order of
   magnitude. A stutter-touching lever runs on it to confirm the mechanism does what it claims;
   its size is measured on real reads only.
6. **One lever at a time on the shipped baseline**; interactions are tested deliberately at
   stage end, not smuggled. One interaction is known in advance: a junk-shape change moves the
   junk-strength optimum, so L1 re-sweeps the two jointly.
7. **The range gate.** No default changes on HG002 alone. Before any keep-as-default verdict:
   the 63-accession tomato cohort at ~3 reads a position runs with the lever on and off, and
   the report states what moved — tract records emitted, the share called heterozygous, the
   QUAL distribution, no-call counts. Tomato has no truth set, so this gate catches behavioural
   breakage (a floor or a discount behaving differently at 3 reads than at 30), not accuracy;
   the report says so each time rather than implying more.
8. **Verdicts are one of three**, written at the end of each section:
   - **default** — ships for every run; requires the range gate and the owner's sign-off;
   - **optional** — reachable via the parameters file or a flag, with one sentence on when a
     user should want it (the fitted-slippage shape is the standing example: right for a run,
     wrong as a constant);
   - **discard** — with the measured number that discards it, so it is not re-proposed.
9. **Prose after the run it reports.** A section quoting a number carries the arm label it came
   from; nothing is written from memory of a superseded run.

---

## 3. Stage 0 — aim before firing (measurements only, no caller changes)

**P0. The baseline pair and the verdict dump.** ✅
Fresh `--defaults` runs at 30× and 50× (the stored callsets predate the 0.05 default), scored
and kept as the program's fixed baseline; the 30× run checked against the sweep's
`outlier0.05` arm, which is the same setting reached the other way. Add a per-tract verdict
output to the instrument (`--verdicts-out`: tract, verdict, error class) so any two arms
crosstab with a join — rule 3 needs it on every arm and today it lives in scratch.

**P1. What the spurious-allele reads have in common.** ☐
At the ~225 spurious-heterozygote tracts, pull the reads carrying the spurious length and ask:
do they share a strand? a start position? a duplicate family? and does the same length persist
in the 300× alignment at the same share? **This is the measurement that aims the whole
program**: clustered reads point at L2 (they are not independent evidence), an even spread that
persists at 300× points at L3/L4 (the locus really yields that length — its own slippage, or
the realigner making it), and no persistence at 300× points at sampling noise no lever fixes.
Deliverable: the 225 partitioned across those three, with cases printed whole for the owner.

**P2. Re-derive the missing-sequence follow-through.** ☐
The split of the never-offered class (434 → aligner / support bar / top-ploidy / unobservable)
predates the corrected comparison and its script has two recorded defects
(`tract_genotype_accuracy_2026-09-03.md` §3.5). Rebuild it on the corrected comparison. This
sizes L4 and gates L7 — the realigner's prize and discovery's prize are both this number.

> **Checkpoint 0:** the three report sections exist; the lever order below is confirmed or
> re-ranked from P1/P2. The owner is told what changed, and rules on nothing unless the order
> moves materially.

---

## 4. The levers, in order

Each is one report section. Build states were verified against the code on 2026-09-03.

### L1 — `junk-shape`: where the junk mass sits *(the owner's pick to go first)* ☐

**Mechanism.** Every read's emission carries a junk term `λ·U`, and `U` is uniform over the
locus's reachable lengths — about 22 of them at a homopolymer, 39 at a five-candidate
dinucleotide (`likelihood/ssr.rs:693`, `ssr_emission.rs:191`). So the floor a read receives is
`λ/22` and `λ/39`: **smallest where the measured junk rate is highest** (1 read in 209 at
period 2+ against 1 in 2,300 at homopolymers), which is why the strength sweep and the literal
reading of λ disagree. Real junk — a chimera, a read from a paralogous tract, a mis-realigned
read — does not land evenly across 39 lengths.

**Arms.** (a) `U` decaying with distance from the nearest candidate, one decay rate, swept;
(b) λ split per period class, two numbers instead of one, set from the measured junk rates and
then swept around them; (c) both. Each arm re-sweeps λ jointly (rule 6). Simulator sanity run:
planted junk reads at a known rate must be absorbed without moving clean loci.

**Target and ceiling.** The thin-evidence tail of the 225 plus part of the 86 "wrong some other
way" — pre-registered exactly in the report before building, but the 158 well-supported cases
are explicitly *not* claimed: a junk term absorbs stray reads, and those reads are not stray.

**Owner input wanted (geneticist's criteria, not a decision):** what junk at a tract should
look like — how far from the true allele a chimera or a paralogous read plausibly lands, and
whether a shape peaked at ±1 unit risks eating real one-unit heterozygotes at low depth.

### L2 — `read-independence`: n agreeing reads are not n pieces of evidence ☐

**Mechanism.** Reads multiply linearly and undiscounted: n identical observations contribute
exactly n times the per-read log-likelihood (`ssr.rs:678,730`), and nothing on the tract path
bounds or discounts anything. If the reads carrying a spurious length share an origin — a PCR
duplicate family, one strand, one mis-realignment — they are one observation counted many
times, and this is the only lever family that can reach the 158 well-supported cases.

**Arms**, cheapest first, the later ones only if P1 says the reads cluster:
(a) an identical-observation discount, `n → 1 + (n−1)·d`, one line at the accumulation site,
`d` swept — this subsumes the owner's `freebayes-cap`, whose `-D` factor is the same idea
applied to the aggregate of unexplained reads (worth at most a tenth, and never bounding one
read); (b) the freebayes form itself, discounting the unexplained mass per genotype;
(c) a per-read floor on how low one read's emission may drive a genotype — the owner's
`gatk-cap`, measured because it costs nothing once (a) exists, with its ceiling pre-registered
as small (it bounds one read; these are many); (d) a beta-binomial overdispersion on the
per-allele read counts — the principled form, and the shape the project already uses on the
SNP path — built only if (a)/(b) show the direction works and leave accuracy on the table.

**Target and ceiling.** The 158 well-supported spurious heterozygotes, *if* P1 shows
clustering; without clustering this lever's premise is false and the section records that and
stops at arm (a).

### L3 — `stutter-em`: the locus's own slippage *(a build: the hook exists, the body does not)* ☐

**Mechanism.** 138 of 242 spurious alleles sit exactly one unit away, concentrated at long
tracts — the signature of a locus slipping far above its stratum's fitted rate. The calling
loop's re-fit is designed (`calling_em_loop.md` §5.1), its configuration and pull-back
constants exist (`SlippageRefitConfig`, 50 pseudo-counts, 20 slipped reads,
`inference/mod.rs:111–280`), and the loop **refuses any non-zero setting** with
`SlippageRefitNotBuilt` — so this is implementing the body inside an interface that is already
there, not plumbing.

**Arms.** The three pull-back settings the spec designs, swept; scored overall **and on the
242-tract subset directly**, because that subset is the hypothesis — if the per-stratum fit
already took what there was, the class will not move and the section says so.
**Interaction to test at stage end:** L3 and the fitted per-stratum rows (L6) overlap by
construction; measure whether the re-fit still earns its keep on top of fitted strata.

**Target and ceiling.** The 138 one-unit cases, minus whatever P1 assigns to the realigner.

### L4 — `realigner`: what the reads are taken to say, upstream of everything ☐

**Mechanism.** Three tracts are verified by hand where the reads plainly carry an allele and
ng's table holds a corrupted spelling of it — a leading base dropped, an interruption
discarded (`tract_genotype_accuracy_2026-09-03.md` §6.2, `src/ng/alignment/ssr_best_path_*`).
This is the only lever that reaches into the 463 never-offered errors, and P1 may show it also
*manufactures* spurious second lengths. It is a defect hunt, not a parameter.

**Plan.** Reproduce `chr3:33,877,690` in isolation — the reads through the tract locus
generator, watching where the sequence is lost; fix; re-run P2's split and the baseline. Sized
by P2 before starting: if the corrected follow-through shrinks the class, this section shrinks
with it.

### L5 — `stutter-keying`: what the stutter model is conditioned on ☐

**Mechanism.** The model is keyed by (period, reference repeat count, slippage group) and
looked up by the candidate's count — candidate-length keying already exists. What it ignores:
the tract's **purity** (interruptions), and the long-tract end where one stratum spans loci
whose rates differ most (the 73-in-1,000 band). Runs only if L3 leaves the long-tract gradient
standing, because a per-locus re-fit and a finer stratification are two answers to the same
observation and the cheaper-to-run one goes first.

**Owner input wanted:** what stratification a geneticist would bet on — purity, motif
composition (A/T against G/C homopolymers), or length alone.

### L6 — `fit-mode`: fitted stutter as something a run can ask for ☐

**Mechanism.** The per-stratum fit from a sample's own reads is worth +0.18 points and cannot
ship as a default — it is a fact about one sample and one chemistry. The machinery to *read* a
fitted parameters file back is built and measured; what is missing is the command that produces
one (`calling_loop_ssr.md` §3.4, deferred). Pure engineering, no research question, and the
program's first natural **optional** verdict.

### L7 — `new-alleles`: discovery, last and gated ☐

**Mechanism.** Admit tract lengths hiding under stutter. The decision half is built and tested
(`calling/inference/discovery.rs`, 12 tests); the wiring belongs inside `select_ssr`
(`tract_genotype_accuracy_2026-09-03.md` §6.5). It targets the **60 collapsed heterozygotes**
and can only enlarge the 225 spurious ones — which is why it runs last: every point L1–L3 take
off the spurious class is headroom this lever gets back.

**Gate.** Runs only if P2's re-derived split still shows a class discovery can reach (the old
count was 61 of 434, superseded) *and* the spurious class has come down. Otherwise discarded
with those two numbers, and the owner — whose lever this is — rules on the discard (§7).

> **Checkpoint 1 (after L1–L3):** compose the winners, re-sweep the traded settings once on
> top, run the range gate, and put the default/optional/discard slate to the owner in one
> sitting rather than lever by lever.

---

## 5. Already answered — not to be re-proposed

Measured under the corrected comparison; report sections exist for each.

| lever | verdict | the number |
|---|---|---|
| `junk-strength` (outlier weight value) | **default 0.05, kept** | thirty-fold range trades 33 spurious for 19 collapsed at homopolymers; 0.05–0.30 flat within 0.1 points; owes the L1 joint re-sweep and the tomato gate |
| flat slip share | **discard** | twenty-fold range, half a point, error classes swing 2× and 9× in opposite directions; shipped 0.10 beside the peak |
| `shorter_share`, `fall_off` | **discard** | every setting within 0.1 points of shipped |
| hand-set length-rising slippage | **discard** | +0.04 at best |
| prior length spectrum (fitting it) | **discard** | the class it would re-balance splits 110 hom-reference against 132 hom-non-reference — a wash by construction |
| prior concentration | **default 1.0, kept** | eighty-fold sweep; 1.0 is the homopolymer peak; +0.04 at period 2+ costs 19 spurious heterozygotes |
| base quality on the tract path | **discard (owner, 2026-09-03)** | informative (35.3 inside against 34.0 outside) but the owner rules the architectural cost out |
| candidate bar, GQ floor, allele-balance rule | **discard** | each costs more than it buys (handoff §5.3), on the old comparison's counts — re-opened only if a lever changes the landscape |

---

## 6. What done looks like

- Every lever above has a report section ending in one of the three verdicts, each verdict
  carrying its number.
- The shipped defaults after Checkpoint 1 are whatever the slate the owner approved says, and
  `defaults.rs` documents each with the run that set it.
- The baseline table in §1 is re-stated at the end with the final defaults, at 30× and 50×,
  both period classes, plus the tomato behavioural comparison.
- Anything discarded is discarded loudly enough that the next session does not rebuild it.

## 7. When the owner is consulted

1. **Any keep-as-default verdict** — the whole slate at Checkpoint 1, with the range-gate
   results beside it.
2. **Discarding a lever from the owner's own list** (`junk-strength`, the two caps,
   `stutter-em`, `new-alleles`) — the number that discards it is presented and the owner rules.
3. **Geneticist's criteria**, asked when the lever reaches them: the plausible shape of junk at
   a tract (L1); what stratification biology favours (L5); and how much heterozygote-calling
   conservatism is acceptable at 3× on a cohort, where every one of these levers pushes the
   same direction — that one question underlies the whole range gate.
4. **Checkpoints 0 and 1**, as marked.

Everything else — arm choices, sweep grids, implementation shape, when a section is finished —
is the assistant's to decide and the report's to record.
