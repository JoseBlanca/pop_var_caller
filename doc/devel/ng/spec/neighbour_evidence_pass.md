# ng — the neighbour-evidence pass: refining each sample's genotypes with the loci around them

*Spec, 2026-11-01. The owner ruled its four first decisions the same day (§9); an adversarial review
the same day found the defects corrected here and the doubts listed in §11, which are the owner's to
rule. It turns the evaluation in [`../research/neighbour_evidence_hmm.md`](../research/neighbour_evidence_hmm.md)
into a design, under the owner's ruling recorded in its §11. It depends on `PL` and `GP` in every
record ([`vcf_output.md`](vcf_output.md) §7.2, branch `pl-gp`) and is to be measured on the founder
cohort of [`str_stutter_benchmark.md`](str_stutter_benchmark.md) §12 (branch `cohort-founders`).*

---

## 1. What this is

pvc calls each locus on its own: a sample's genotype comes from its reads at that locus and a prior
from the cohort's allele frequencies. In a cohort of inbred lines the loci around it say more. The
lines descend from a few haplotypes, so a sample's chromosome is, stretch by stretch, a copy of
other samples'; and an inbred sample's two copies are one copy counted twice over long runs. The
pass uses both: for each sample, it asks along the chromosome *which other samples is this one
copying here, and are its two copies the same here?*, and re-calls each genotype from its own reads
together with the answer.

**What the prototype measured** (report §4–§5), on tomato1's 57 accessions over 80 regions of
100 kb, against the nine deep accessions' own full-depth calls:

- with the nine cut to 3× and 1× among 48 accessions at native depth, 46% and 56% fewer genotype
  errors; with the whole cohort cut to 3× and 1×, 48% and 36%;
- at sites where the shallow accession had no read, 85% of true non-reference genotypes wrong today
  and 26% with the model — when the nine sat among native-depth accessions; with the whole cohort at
  1× the same comparison is 88% and 47%;
- in a simulated outbred cohort, nothing.

Those numbers were measured with a transition between the autozygous states that this review found
wrong (§3.5) and with a procedure that differs from this spec's in the ways §3.7 lists; the Rust
pass must re-measure them (§8).

### 1.1 Goals

- Re-call every sample's genotype at every PASS record pvc writes, SNPs, indels, multi-allelic
  records and STR records alike, from the sample's own likelihoods and the cohort's shared haplotypes.
- Give a genotype, with an honest quality, where a sample has no reads but its neighbouring sites
  have its evidence (§4.2).
- Run on by default where it helps and off where it cannot, saying which and why.
- Work from what a record carries — `PL`, `GP`, the alleles — and each sample's inbreeding
  coefficient from the parameters file, so the same code can run on a finished VCF (§9 Q4).

### 1.2 Non-goals

- **Phasing — in two companion specs.** The pass phases samples to give each other haplotype
  templates ([`neighbour_evidence_phasing.md`](neighbour_evidence_phasing.md)) and writes the phase
  in confident blocks ([`neighbour_evidence_phase_output.md`](neighbour_evidence_phase_output.md))
  (owner, 2026-11-02).
- **New loci.** A locus pvc did not build stays unbuilt; the pass changes genotypes only. On one
  tomato accession alone at 3×, 42% of the true non-reference genotypes were at such loci (report
  §5.2); that loss belongs to candidate generation, not here.
- **Read linkage** — reads covering two sites. Left out by the owner's ruling (report §11).
- **Outbred cohorts, without phased templates.** Using shared haplotypes there needs phased templates
  and repeated sweeps over the cohort
  ([`07_decision`](../research/phasing_imputation/07_decision_separate_tool_and_route.md) §2); this
  spec's pass does neither and is off there (§5). Partly outbred cohorts are what
  [`neighbour_evidence_phasing.md`](neighbour_evidence_phasing.md) adds.

### 1.3 Vocabulary

- **Template.** Another sample, as a source a copy can copy from. At each record a template is its
  *allele shares*: the expected fraction of its copies carrying each allele, from its `GP`. An inbred
  sample's shares are near 0 or 1, which is why a genotype can stand in for a haplotype.
- **Autozygous.** Of a stretch where a sample's two copies are one ancestral copy counted twice; the
  sample is homozygous along it.
- **Informative record** of a sample: one where its `PL` is not flat — its own reads said something.
- **Window.** A run of consecutive records the pass holds and solves together (§2.2).

---

## 2. Where it sits

### 2.1 A stage over finished records

**The pass reads finished VCF lines and writes them rewritten.** With the hidden-duplication filter
on, records are parked on its spill as finished lines and read back by its third pass; with it off,
they are encoded for the writer, and the stage takes them as the lines the writer would write. One
form for both runs and for the stand-alone command (§9 Q4) means one code path, and the same input
precision everywhere: `PL` in whole Phred units, `GP` at three decimals.

This is the "pass over the emitted records" [`run_streaming.md`](run_streaming.md) §4.3 requires for
any coupling between segments: nothing about calling changes. **The stage's input is a bounded
buffer**: if the pass falls behind, the calling workers wait for it rather than the buffer growing
without limit; how often that happens is measured (§8).

**Built (2026-11-01) as the tightest such buffer: the stage lives inside the VCF writer** and runs
on the thread that writes. Every line, with the filter on or off, already ends in one call, the
writer's `write_line`, made on the calling thread; the writer checks the line's place against the
one before it as before, hands the line to the stage, and writes whatever the stage hands back;
finishing the file finishes the last contig's windows. So calling waits while a window is solved,
the stage never holds more than a window, and neither the filter's third pass nor the filter-off
path changed. A thread of its own, letting calling run ahead by a bounded number of lines, is the
change to make if §8.6 finds the wait matters (`src/neighbour_evidence/stage.rs`).

### 2.2 Windows

The pass cannot write a record until it has seen the records after it, because a sample's state at a
site depends on both sides. It holds records in **windows**: a *core*, written when the window is
solved, with *margins* either side read for context and written by the neighbouring windows.

- **Which record belongs to which core is decided by arrival order**, not by position: each window's
  core is a run of consecutive records, so every record is written exactly once and in the order it
  arrived.
- **Positions are each record's written `POS`**, asserted non-decreasing within a contig; distances
  are differences of `POS`.
- **Windows span the gaps between BED regions.** A gap is only a long distance, and the corrected
  transition (§3.5) returns the state to its long-run share across it. A window is cut only at a
  contig's end; nothing crosses contigs.
- **Sizes are set by length and by a memory budget in bytes** (§6): the core starts at 2 Mb and each
  margin at 1 Mb, and both shrink as the cohort grows, the margin never below 500 kb. The core
  shrinks first, down to the margins' length; then core and margins shrink together, in steps of
  50 kb, to the floor. A window's size is chosen from the bytes its records would hold, so a window
  waits until the records its largest affordable size reaches have arrived
  (`src/neighbour_evidence/windows.rs`).
- **A window over budget even at the floor stops the run** with an error naming the window, its
  records, samples and bytes, and the budget to ask for (`--refine-window-bytes`). The run neither
  switches refinement off nor goes on past the budget (owner, 2026-11-13): a run that cannot refine
  says so, rather than writing a file whose header promises a pass it did not run, or running out of
  memory ten minutes into a window.

### 2.3 Parallelism

Within a window, every sample is solved on its own once the templates are fixed (§3.3), so the
samples of a window are shared among the run's threads. Windows are solved one after another, which
keeps the written order without reordering. **The output does not depend on the thread count**:
each sample's walk is a function of the window alone, and template ties break by sample order (§3.3).

---

## 3. The model

### 3.1 In words

For one sample, walk along the window's records. At each, the hidden state says either *the two
copies are copying templates `k₁` and `k₂`*, or *the two copies are one, copying template `k`*.
Between neighbouring records a copy keeps its template with a probability that falls with the
distance between them, and the autozygous state persists the same way. At each record, a state is
scored by how well the sample's reads fit the alleles its templates carry. A forward and a backward
walk give, at every record, the probability of each state given the whole window; the genotype
posterior is the sample's own likelihood weighted by what the states say.

### 3.2 The inputs, read off the record

For a record with `A` written alleles and, at diploidy, `G = A(A+1)/2` genotypes:

- **The sample's likelihood** `L(g)` for each genotype, from `PL`: `L(g) ∝ 10^(−PL_g / 10)`.
- **Each template's allele shares** `q(a)`, from its `GP`: `q(a) = Σ_g GP_g · copies_a(g) / 2`.

What a record contributes depends on what it carries:

- **A record whose `FILTER` is not `PASS`** — tagged by the hidden-duplication filter, or `EMNoConv`,
  `tooManyAlleles`, `lowDepth` — **is passed through unchanged and is no site of any walk.** A
  collapsed paralog's false heterozygotes would otherwise break the autozygous runs around them.
- **A record with no rows at all** — one merged from an STR record and its flank's generic records
  (`vcf_output.md` §7.2) — is passed through unchanged and is no site.
- **A sample written `.` in both** (set aside by the candidate step) is not re-called at that record,
  and as a template there contributes equal shares, `1/A` each — no information.

### 3.3 The templates

For each sample and window, the `K = min(16, N − 1)` other samples that most often disagree with it
least **over short stretches of the window, against what the cohort's allele frequencies predict for
an unrelated sample, judged only where both have reads, with more data preferred between two that
disagree about as little** (owner, 2026-11-02 and 2026-11-05). The records compared are the sample's informative records in the
window — a fixed sample of at most 500 of them, spread evenly over them by record order (the `i`-th
of 500 is record `⌊i · n / 500⌋` of the `n`), so the choice is deterministic and weighs the whole
window, not its left margin. Of those, a candidate is compared where it is informative too. At each
such record:

- the **observed** disagreement is the chance that one allele drawn from each sample differs,
  `1 − Σ_a q_self(a) · q_other(a)`;
- the **expected** disagreement is the same against the cohort's frequencies,
  `1 − Σ_a q_self(a) · p(a)` — what an unrelated sample would show.

**Judged stretch by stretch, ranked by stretches won (owner, 2026-11-05).** The compared records
are cut, in record order, into **stretches** that each hold `I` units of expected disagreement: a
stretch ends at the record where `e`, summed over its records, reaches `I`, counting only records at
which some candidate is informative; what is left after the last such record joins the stretch
before it, so no stretch is judged on less, and a window holding less than `I` is one stretch. One
unit is about two SNPs at frequency ½ where the sample is certain, or about twenty at 5%.

**`I` grows with the cohort (owner, 2026-11-06):** `I = max(1, ln(N / 8))`, `N` the other samples —
one unit up to 21 of them (the 22-sample F2), 1.6 at 39, 3.9 at 399, 4.8 at 999, 6.4 at 4,999 (`templates::stretch_information`). `N` counts samples in every round, also from round 2,
where each offers two haplotypes but takes at most one of a stretch's places. In each stretch every candidate gets a score from its own `d`
and `e` there:

- the score is the mean of the rate's Beta(`d + 1`, `e − d + 1`) posterior plus two of its standard
  deviations — pessimistic, so a rate measured on little scores worse — with `d` capped at `e`, so a
  candidate less alike than an unrelated sample scores as one;
- a candidate with `e` = 0 there — no shared record, or only records where every sample carries the
  same allele — has nothing to be judged on and cannot win the stretch;
- the two lowest scores **win** the stretch, and so does any candidate tied with the second. From
  phasing round 2, where the candidates are haplotypes, the two places go to two samples, each
  scoring by its better copy — so a homozygous sample's two equal copies do not take both
  (phasing spec §2.3).

Candidates rank by stretches won, most first; then by the same score over the whole window's
compared records; then by sample order. The first `K` are the templates. **A template's stretches
won also set where a copy that changes template lands** (§3.5): with weight `(wins + ½)⁴`, so a
template that won twice as many stretches is landed on about 16 times as often, and one that won nothing
is still reachable, at an 81st of the weight of one that won a single stretch.

**Why stretches — changed 2026-11-05.** A sample's two copies change template along the window at
every crossover in its ancestry. A whole-window score ranks highest whoever is fairly close
everywhere; the samples a copy actually copies are close over a part each — an F2 plant's parents
over the stretches either side of a crossover — and only stretches short enough to fall between
crossovers see that. On the simulated F2 (two inbred parents at 30×, 20 offspring at 1×) the
whole-window rule often left an offspring's parents out of its templates or buried them among
siblings; stretch wins raised the heterozygous offspring genotypes called right from 72 and 67 in
100 to 82 and 81 on two simulations, with the wrong calls among them falling from 6 and 7 in 100 to
3. In 40-sample natural cohorts at 1×, heterozygotes right rose by about 5 in 100 (outbred) and 4 to
6 (selfing); samples at 5× gained up to 3 in 100, and samples at 10× and above did not move
(`doc/devel/reports/implementations/local_sources_2026-11-05.md`).

- *Why two winners:* tuned on those cohorts; four measured worse, and so did fixed stretches of 50
  records and growing a stretch until its candidates separated. Two winners are one a copy.
- *Why `I` grows with the cohort — changed 2026-11-06:* among many candidates some unrelated one
  matches a short stretch as well as a relative by chance, and takes the relative's place. At 400
  samples one unit gave a stretch 6.5 winning samples where two were asked for, against 2.45 at 40,
  and the phase lost about a tenth: heterozygotes phased right fell from 108,656 to 98,871 of
  145,951 in the refine_sim benchmark's outbred cohort. Longer stretches there recovered it and
  called fewer true variants wrong — 4 units: 108,600 phased right and 2,201 wrong against 2,771; 6
  units: 111,110 and 2,145. The F2 runs the other way: on 20 offspring at 1×, every length above one
  unit called more heterozygotes wrong or missing (266 wrong with reads at one unit, 298 at two, 337
  at three), and making its crossovers ten times rarer did not change that, so it is not the
  benchmark's dense crossovers; why the F2 wants short stretches is not traced. At 40 samples
  (outbred) and 400 (selfing) lengths 1 to 6 measured level, within 1 in 100 of every genotype class.
  The form follows a crude count — if an unrelated candidate matches `I` units with chance about
  `e^−I`, the cohort's chance matches are `N·e^−I`, and `I = ln(N / 8)` holds them to about eight —
  and the 8 is measured, the value that keeps the 22-sample F2 at one unit, its best length (3,087
  heterozygotes phased right there, 2,837 to 3,043 at the other lengths from 0.8 to 1.5), and gives
  400 samples 3.9
  (`doc/devel/reports/implementations/stretch_scale_2026-11-06.md`).
- *Why land by wins:* on the F2, landing every template alike kept 74 heterozygotes in 100 right
  against 81; in the natural cohorts it made no difference.
- *Why to the fourth power — changed 2026-11-06:* a 1× F2 plant's parents rank first, but its
  siblings still win about two of every three stretch places: at one read a site, a heterozygous
  stretch looks as much like a heterozygous sibling as like either parent, and no single candidate
  can show it is the pair of parents. Landed on in proportion to their wins, the siblings took
  copies the parents explain with certainty, and their soft shares left heterozygotes uncalled.
  With the parents as the only templates the local F2 called 89.7 heterozygotes in 100 right
  against 85.2; landing at the fourth power gave 89.5, without being told who the parents are. On
  the refine_sim benchmark's F2 it phased 4,607 heterozygotes right against 2,775, missed 208
  variant genotypes against 284 and left 795 genotypes without reads uncalled against 1,050, for 19
  more false variants (13 → 32); the 400-sample outbred and selfing cohorts measured level, within 1
  in 100 of every class (outbred: 61 fewer missed, 63 more false). Powers 4 to 6 measured alike on
  the F2; 2 gained nothing (`doc/devel/reports/implementations/landing_power_2026-11-06.md`).
- *Tried and not adopted* (same report): IMPUTE5's choice, the haplotypes sharing the longest exact
  run of likeliest alleles with the sample, in rounds 2 and later — no better than stretch wins on
  the natural cohorts, worse for selfing at 1×, and adding nothing measurable on top of them; and
  Beagle's composite templates, one re-chosen in every stretch — re-chosen from the sample's own
  reads every few records, they echo those reads back, and the F2 fell to 18 in 100.

**Data wins near ties, and certainty is not penalised.** Within a stretch the pessimistic score
charges a candidate for how little it is judged on, so a deep sample beats a shallow one that agrees
as well — and it is the deep sample whose own data fills a site the target has none at. **The cost,
stated:** a shallow candidate cannot win a stretch it has no reads in, so a shallow close relative
wins fewer stretches than a deep sample that is only fairly close, and may rank below it. Where that
loses accuracy — a cohort with shallow close relatives and many deep distant samples — is not
measured.

**Above 400 samples — not measured.** `I` was set from cohorts of 22, 40 and 400 samples. At
several thousand it gives 6 to 7 units, about 12 to 14 SNPs at frequency ½; whether that keeps
chance matches out without stretches spanning a copy's change of source is not measured. A cohort
of close relatives that is also large — a big family or breeding population — gets long stretches
though the F2 measured best at short ones; that case is not measured either.

**Why against the frequencies — changed 2026-11-05.** The rule before scored the plain
disagreement, `½ Σ_a |q_self(a) − q_other(a)|`, over the shared records, and it favoured exactly
the samples that carry least information:

- a sample whose genotype is uncertain has shares near the frequencies, and two uncertain samples
  differ little from each other — so a 1× sibling looked closer to a 1× F2 plant than its 30×
  parents, whose shares are 0 or 1;
- a record where the whole cohort carries one allele is agreement everyone shares, and counted as
  much as one that separates relatives. A 1× sample's confident records are mostly of that kind.

On the simulated F2 (two inbred parents at 30×, 20 offspring at 1×), most offspring copied neither
parent, and 30 of 100 heterozygous offspring genotypes came out right; with this rule and the
background as a last resort (below), 51 of 100. In two 40-sample natural cohorts with samples at 1×
and 2×, both changes together gained 3 to 7 genotypes in 100 at 1× and did not move the samples
at 5× and above by more than 1 in 100
(`doc/devel/reports/implementations/refine_sources_2026-11-05.md`).

Why only where both have reads: a sample without reads carries the cohort's frequency, and counting
that as disagreement would rank a shallow close relative below a deep distant one. On the founder
simulator the change from the earlier rule, which did count it, moved fills and errors by at most a
dozen genotypes in any cohort, uniform or mixed in depth (10 samples at 30× among 40 at 3× or 5×; 5
among 15), because a 3× sample has reads at about 95% of records and the earlier rule's charge for
the rest was small. It is kept as the rule that says what it measures; where it should matter — 1×
shallow samples among deep ones, cohorts far larger than `K` — is not measured. **The sample itself is never a template.** `K` = 16 was
what the prototype ran; 2, 4 and 8 measured worse on every cohort of the 2026-11-05 report.

**A background template** joins the `K` (§11 Q6): its shares at each record are the cohort's allele
frequencies — `AF`, the reference taking the rest. A copy that matches none of the sampled
haplotypes — an introgression, a private lineage — copies the background, and its genotype falls back
to what the cohort's frequencies say rather than to the nearest template's allele. It takes part as
a template does — pair states with it, an autozygous state copying it, the same mismatch blur — but
**it is a last resort: a copy that changes template lands on it with probability 10⁻⁵** (§3.5;
owner, 2026-11-05). As one template among 17 alike it was copied in place of real ones wherever an
allele stood near 50%: at one read a site, a source that predicts each allele with even odds explains
a heterozygote as well as the two haplotypes that make it, and in an F2 every allele stands at 50%.
At 10⁻⁵ it is copied only where no sampled template explains the reads; across the cohorts measured,
any share at or below 10⁻⁴ gave the same results as leaving it out, so the share's exact value does
not matter there. The stretches it exists for — a haplotype no other sample carries — are absent
from those cohorts, so how well it still serves them is not measured.

A sample with no
informative record in the window has no basis for choosing templates and is not re-called in it.

**A known circularity, stated.** A template's `GP` came from a prior built from the other samples'
allele copies — the target's included. So a little of the target's own evidence returns through its
templates. Its size is not measured; it shrinks as the cohort grows, and §8 measures it at 10 samples
by re-calling with templates whose `GP` were computed without the target.

### 3.4 The emission

A template's shares are blurred by a mismatch rate `ε` — a mutation since the shared ancestor, or the
template's own wrong call — `q̃(a) = (1 − ε)·q(a) + ε / A`, with `ε` = 0.01 as in the prototype. Then:

    not autozygous, templates (k₁, k₂):   e = Σ_a Σ_b q̃_k₁(a) · q̃_k₂(b) · L(a/b)
    autozygous, template k:               e = Σ_a q̃_k(a) · L(a/a)

The sums run over alleles, so a record with eight alleles, an STR among them, costs more arithmetic
and nothing else.

### 3.5 The transitions

Between records `d` bases apart:

- **A copy keeps its template** with probability `e^(−d/S)`, and otherwise draws one: the
  background with probability 10⁻⁵, each of the `K` sampled templates a share of the rest in
  proportion to its stretches won plus ½ (§3.3; before 2026-11-05 every one of the `K + 1` alike). In the non-autozygous state the two
  copies do this independently.
- **The autozygous state is a two-state Markov chain** with runs of mean `R` and gaps of mean
  `R′ = R · (1 − F) / F`, so a share `F` of the sample is autozygous in the long run:

      stay autozygous   = F + (1 − F) · e^(−d/τ)
      stay not          = (1 − F) + F · e^(−d/τ),        1/τ = 1/R + 1/R′

  Across a gap far longer than `R` and `R′` — a centromere, the space between BED regions — the state
  returns to `F`. *Corrected 2026-11-01:* the first version, and the prototype the report measured,
  used `e^(−d/R)` and `e^(−d/R′)` as the two stays, which across a long gap swaps the states: a
  sample 99% autozygous before a centromere came out 1% after it.
- A change of state draws its templates afresh, with the same landing chances, and the walk's first
  record starts from them too: autozygous on template `k` with `F · w(k)`, the pair `(k₁, k₂)` with
  `(1 − F) · w(k₁) · w(k₂)`, `w` the landing chance (`WalkParameters::landing_weights`,
  `src/neighbour_evidence/model.rs`; the weights from `ChosenTemplate::landing_weight`,
  `templates.rs`).

**The parameters, per sample.** `F` is the sample's fitted inbreeding coefficient from the
parameters file, clamped to [0.01, 0.99], and the walk keeps it (changed 2026-11-06). It was
re-estimated as the mean autozygous posterior over the sites the fitting windows walked, and that
ran away where every sample is shallow: at 2×, one soft template explains both of a sample's copies
about as well as two do, and in a 40-sample outbred cohort all at 2× the re-estimate reached 0.83 on
average, up to 0.99, against a run-level estimate of 0 to 0.4 and a truth of 0. The pass then
called heterozygotes homozygous — 3,344 true heterozygotes right against 7,521 with refinement off;
keeping the run-level `F` gave 5,007. On every other cohort measured (inbreeding 0, 0.5, 0.8 and
0.9, 2× to 30×, and the F2) keeping it changed no genotype class by more than 10 in 17,000, and
even a selfed cohort (F 0.9) walked at `F` 0.01 moved by at most 7: the walk sees autozygosity
through its templates (`doc/devel/reports/implementations/fixed_f_2026-11-06.md`).

`S` and `R` are chosen by likelihood from small grids — `S` ∈ {30, 300} kb,
`R` ∈ {10, 30, 100, 300, 1,000, 3,000} kb — **on the first windows of the run that hold at least
5,000 of the sample's informative records, and then frozen** for the rest of the run. Until then
each window adds its twelve walks' log-likelihoods to a running total a grid point, walks with the
point whose total is largest (ties to the smaller `S`, then the smaller `R`) — so a second fitting
window refines the first one's choice rather than replacing it. A fitting window costs thirteen walks: the grid's twelve and the one that gives its posteriors.
*Settled in code 2026-11-01 (`src/neighbour_evidence/fit.rs`); the sentence before left open how
several windows combine.* Refitting every
window would multiply the pass's cost by the grid's twelve walks (§6); whether a frozen fit loses
accuracy is a measurement (§8).

### 3.6 The posterior

At a record, with `P(state | window)` from the two walks and `e(state)` the record's emission:

    P(g | window) ∝ L(g) · Σ_states [ P(state | window) / e(state) ] · P(g | state)

where `P(a/b | k₁, k₂) = q̃_k₁(a) q̃_k₂(b) + q̃_k₁(b) q̃_k₂(a)` for `a ≠ b` (and the one term for `a = b`),
and `P(a/a | autozygous k) = q̃_k(a)`. Dividing by the emission removes the record's own reads from the
state probability, so they enter once, through `L(g)`. The sum over pair states is computed as
`Q̃ᵀ · W · Q̃` — templates by alleles, then alleles by alleles — at `K²·A + K·A²` operations a record,
not `K²·A²`.

### 3.7 Where this departs from the prototype

The numbers in §1 were measured with a procedure that differs from this one, so a gap between the
Rust pass and those numbers is not by itself a defect:

| | prototype (`eval_tomato.py`) | this spec |
|---|---|---|
| templates chosen | once a sample, over all its records in the 8 Mb, flat ones included | per window, over its informative records |
| `R` and `F` | from a separate two-state model with calling's prior branches | inside the copying model |
| `S` | chosen per sample from {30, 300} kb | the same, on the first windows, frozen |
| the autozygous transition | the swapping form (§3.5) | corrected |
| set-aside samples | re-called, with a prior from the others' shares | not re-called |
| non-PASS records | sites like any other | passed through, no sites |

---

## 4. What it writes

### 4.1 The fields

| field | after the pass |
|---|---|
| `GT` | the argmax of the new posterior, or `./.` where §4.2 withholds it |
| `GQ` | from the new posterior, as [`calling_quality.md`](calling_quality.md) §4 defines it |
| `GP` | the new posterior |
| `PL` | unchanged — it is the evidence |
| `DP`, `AD` | unchanged |
| `REPCN`, `REPCNB` | each called allele's value, read off another sample's column that names the same allele |
| `INFO` `AC`, `AN` | recounted from the new genotypes |
| `QUAL`, `AF`, every other `INFO` field, `FILTER` | unchanged |

- **`REPCN` and `REPCNB` need no tract arithmetic.** Calling writes only alternatives some sample
  carries (`drop_alternatives_no_sample_calls`), so every allele of an input record is named by some
  input call, and its repeat count is in that sample's column.
- **Alleles the pass leaves uncarried are dropped** as calling drops them: the alternative goes from
  `ALT`, and every sample's `PL` and `GP` keep only the genotypes over the written alleles
  (`vcf_output.md` §7.2's subsetting). **A record left with no alternative carried is dropped** — the
  same rule calling applies to a locus whose alternatives no sample is called with. *Both are the
  owner's to confirm, §11 Q2.* As in calling, a record where no sample is called at all keeps its
  alleles. The subsetting reaches every per-allele field — `AD`, and `AF`, `MQALT` and `MQDIFF` in
  `INFO` — and re-shifts `PL` so the best kept genotype is 0. **One thing it cannot redo**: calling
  drops alleles before it widens a repeat tract's alleles over their flanks, so a tract record that
  loses an allele here keeps the padding computed over the larger set — a longer `REF` than calling
  would have written, still a correct description of the alleles.
- **`GQ` comes from the written `GP`** — after the subsetting, by the same function calling uses
  (`calling_quality.md` §4); a withheld sample's `GQ` is `.`, as calling writes for a no-call.
- **A repeat count no input column names is written `.`.** The rule above holds for alternatives;
  the reference allele can be named by no input call (every sample homozygous for an alternative),
  and re-calling a sample to carry it then has no column to read its count from. The run's log
  counts these.
- **The header** gains one line recording whether the pass ran and why, and its parameters, and `GP`'s
  declaration says the posteriors are the pass's. A run where the pass is off writes the line too,
  with `ran=no` and the reason, and no parameters:

      ##neighbourEvidence=<ran=yes,why="57 samples, median inbreeding 0.98",templates=16,mismatch=0.01,core=2000000,margin=1000000>
      ##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype posterior probabilities after the neighbour-evidence pass, VCF genotype order">

### 4.2 A sample whose reads said nothing

A sample with flat `PL` at a record — written `./.` by calling (`vcf_output.md` §7.1) — gets the
pass's posterior, and **is written with a genotype when the pass gives one from evidence**: its own
informative records at the neighbouring sites. Its `GQ` says how sure that is, and its flat `PL` says
none of it came from its own reads (§9 Q1).

**But a posterior with no evidence behind it is not written as a call.** Where a sample has no
informative record nearby, every emission around the site is 1, the state posterior is its prior,
and the posterior genotype is the cohort's frequency — `0/0` at `GQ` 17 for a 2% allele, the
confident prior-only call `vcf_output.md` §7.1 exists to forbid. So **a flat-`PL` sample is written
with a genotype only if it has an informative record within one segment length `S` on each side of
the site** (on one side, where the contig ends within `S` on the other — at a contig's right end only
when the header gives its length); a record at the same position counts on neither side, since it
shares the site's reads rather than bracketing it; otherwise it stays `./.` and
its `GP` is the pass's posterior, as calling writes a prior. **And only if the call's `GQ` is at
least 9** (owner, 2026-11-02 at 10; 9 since 2026-11-08): below it the sample stays `./.` with the
pass's posterior as `GP`. The floor is set for cohorts of 100 samples and more, where the
refinement is meant to run (owner, 2026-11-08). There, calls without reads at `GQ` 9 were wrong 2
to 6 times in 100 in outbred cohorts and an F2 of 100 plants, and 9 in 100 in a selfing
population. 8 was not taken, the conservative choice: its calls were wrong 3 to 8 times in 100 in
the outbred cohorts and the F2, but 31 in 100 in the selfing population. Lowering the floor from 10
to 8 cut the missing SNP genotypes of the shallowest samples by a fifth to a third for 1.5% to 4%
more SNP errors
(`doc/devel/reports/implementations/no_read_fill_gq_2026-11-08.md`). In smaller cohorts the floor
was set at 10 because, on the founder simulator at 1× (30 samples), calls below it were wrong 28
times in 100, mostly a common alternative called homozygous in a sample carrying the reference;
evidence nearby that does not say which samples a sample copies leaves its posterior near the
cohort's frequency (`doc/devel/reports/implementations/neighbour_pass_phase4_2026-11-02.md`,
decision 1). *This narrows ruling Q1 to the cases its
reason covers; the owner's to confirm, §11 Q1.*

**And only if the sample is copying the background template there with probability at most 0.1**
(owner, 2026-11-02): the walk's posterior probability that the autozygous state copies the
background, or either copy of a pair does (§3.3). Where no other sample carries the sample's
haplotype, the call is the cohort's frequency again, however confident its `GQ`. On the founder
simulator's inbred cohorts, the calls this withholds between 0.1 and 0.5 were wrong 13 times in 165,
and those above 0.5 85 times in 719, against 1 to 2 in 100 for the calls kept; at 10 samples and 1×
the limit keeps 101 of 294 calls, 2 of them wrong, from 47 (same report, the background-copy
section).

### 4.3 Each sample's fills, checked on its own data

The test of hiding genotypes the reads settled and asking the pass for them is done for every
sample at every such record, at no cost (owner, 2026-11-02): the walk's prediction without the
sample's own reads at a record is what a no-read call there would be written with. At each core
record where the sample's own evidence settles its genotype — posterior at least 0.999 from its
reads and calling's prior, the cohort's allele frequencies and its inbreeding — and where the §4.2
rules would have written the prediction as a call, the prediction's genotype is compared with the
settled one. The prior has to be in "settled": at 3× reads alone never tell a homozygote from a
heterozygote with confidence, and a check that needed them ran on no shallow sample — the samples
that are filled.

The run reports, for each sample, the no-read calls written, the records checked and how many the
prediction got wrong: in a table beside the VCF, `<name>.fill_check.tsv`, and summarised in the log,
with the share given only for samples checked at 30 records or more. It is a best case — a record
no read reached can be harder than one where the sample had reads.

On the founder simulator (`validate_fill_check.py`): at 3× the cohort's estimated fill error rate
against its true one was 1.10% against 1.52%, 0.63% against 0.49%, 0.13% against 0.15%, and on the
outbred cohorts 4.16% against 2.90% and 5.52% against 4.55%; samples reported at 0–0.5% were truly
wrong at 0.17% of their fills, at 2–5% at 1.5%. At 1× almost no sample has a settled genotype and
the share is not given — the corner where fills are worst (2–2.5% wrong) goes unestimated, and the
log says so.

---

## 5. When it runs

**On by default when the cohort can supply templates and is inbred enough for genotypes to stand in
for haplotypes; off otherwise, saying so in the header and the run report.** A run-level decision,
taken once, before the first record is written:

- **At least `N_min` samples.** Starting value 10. On the simulation, 10 inbred samples from 8
  founders already removed 46% of the errors at 3×, against 72% at 63 samples (report §4); at 10
  samples `K` is 9.
- **The cohort's median fitted inbreeding coefficient at least `F_min`.** Starting value 0.5. The
  outbred simulation gained nothing; between the two ends is unmeasured.

Both are starting values, to be set from §8's measurements. `--refine-genotypes auto|on|off` (named
`--neighbour-evidence` until 2026-11-03),
`auto` the default, overrides the decision either way — except that `on` cannot run the pass on one
sample or at a ploidy other than two. The median of an even count is the mean of the middle two.
Within a window a sample's `K` is fixed at `min(16, N − 1)`; a template without a row at a record
carries equal shares there (§3.2), and a record where fewer than two samples have rows is no site.
The decision goes into the run's log as well as the header, and after the run one log line says
what the pass did: records re-called, genotypes called (and how many of them where the sample had
no reads), withheld, changed, alleles and records dropped. **The run's record count is read from
the file**, since the pass can drop records.

**At a ploidy other than two the pass is off**: the pair state is diploid.

---

## 6. One sample and three thousand, three reads and three hundred

- **One sample.** Off: there are no templates. A lone sample would gain 5% at 3× on tomato from the
  autozygous state alone (report §5.2), and the owner ruled that not a goal.
- **Ten samples.** `K` is 9; templates are few, and the circularity of §3.3 is largest. §8 measures
  both.
- **Three thousand samples.** Own arithmetic at tomato's density (26 records a kilobase in tomato1,
  so about 24 million over a 900-Mb genome), `K` = 16, a biallelic record:
  - *Choosing templates*: `N²` comparisons over 500 records a window, 4.5 × 10⁹ a window; at the
    1.5-Mb windows this cohort size needs (below) about 1,800 windows, 8 × 10¹² over the genome —
    about half an hour.
  - *The walks*: forward and backward cost about `2 · 3 · (K² + K)` ≈ 1,600 operations a sample a
    record, and the posterior `K²·A + K·A²` ≈ 580. Every record is walked in its core and in its
    neighbours' margins: twice at the default sizes (a 2-Mb core, 1-Mb margins), three times at the
    floor this cohort size reaches (0.5-Mb core and margins) — about 6,500 a sample a record. Over
    24 million records and 3,000 samples that is 4.7 × 10¹⁴ operations, **about 33 core-hours** at
    4 × 10⁹ a second; at 57 samples, at default sizes, 25 core-minutes. The parameter grid adds its twelve walks on the first windows only (§3.5).
  - **Measured 2026-11-13, about 150 times that.** On 100 to 800 tomato25 samples over 1 to 4 Mb of
    chromosome 0, 48 threads, the cost grew in proportion to the samples (773, 1,415, 2,871 and
    5,281 CPU-seconds for 100, 200, 400 and 800 over the same 1 Mb), with no sign of the `N²`
    template choice: 3 CPU-seconds a window at 400 samples. The walks were 85% of the time, about
    130 µs a sample a record in a window after the fitting ones. Two factors make up the gap: the
    phasing rounds (`neighbour_evidence_phasing.md`) walk five times, the later four drawing nine
    paths each with reads, about eight walks' work; and a walk ran about 20 times slower than the
    operation count above (12 µs a sample a record measured, 0.55 µs counted). At 2,176 samples that
    projects to about 4,900 core-hours for the genome, five days on 48 cores
    (`doc/devel/reports/implementations/refinement_memory_2026-11-13.md`; making the walk cheaper is the next
    step).
- **Three reads a position.** Where it helps: the prototype's gains are all at 1× and 3×.
- **Three hundred reads.** The sample's own likelihood dominates every emission and the posterior is
  the calling's; the pass costs time and changes nearly nothing. It still runs if the switch says so:
  depth does not turn it off, because a deep cohort can hold shallow samples.

**Memory.** A window's memory has two parts: what it holds for every sample at every record, which
grows with the cohort, and what each thread's walk holds at every record, which grows with the
threads.

- **Per sample and record**, for every record the window spans: the record's VCF line, held
  compressed; its likelihoods and shares (`L` and `q`, as `f32`); each sample's read alleles, seven
  bytes each (`PackedRecordFragments`); and each sample's two haplotypes, the templates they
  followed and its phase across the rounds. For each record of the core, besides: each sample's sums
  over the averaged rounds and what it writes. *Measured 2026-11-13* on 400 tomato25 samples
  (chromosome 0, about 22 read alleles a sample a record there), about 140 bytes a sample a record
  held before the solve and about 100 more in it, biallelic; *before* that day's changes about 500,
  of which the read alleles at twelve bytes each were about half. The spec's first estimate, 40
  bytes, counted only the likelihoods and shares.
- **Per thread and record**: the walk's forward probabilities, `K² + 3K` numbers, its terms for the
  backward sweep, the round's inputs and, from round 2, the nine drawn paths — about 4 to 5 KB a
  record, so 48 threads over a 70,000-record window hold about 15 GB. Its emissions, `K² + K` numbers
  a record, are kept for two records only and computed again in the backward sweep (2026-11-13).

**The window is sized to a budget in bytes**, `--refine-window-bytes`, 32 GB by default. It counts
the per-sample part exactly as the records arrive — each line as it arrived, though most are held
compressed, so the count does not depend on when its batch was compressed. It does not count the
per-thread part: that would make the windows, and so the output, depend on the thread count (§2.3).
At 2,176 samples and about 250 bytes a sample a record, 32 GB holds a window of about 58,000
records, 2.2 Mb at 26 records a kilobase; the 500-kb floor is reached at about 3,000 samples, and
past it the run stops (§2.2).

---

## 7. Traps

- **Genotype order.** `PL` and `GP` are in the shared table's order (`vcf_output.md` §7.2); the pass
  must read and write them through the same table, never by a hand-rolled index.
- **A moved `GT` moves `REPCN`, `REPCNB`, `AC` and `AN`** (§4.1); a stale `REPCN` is a silent error.
- **A flat-`PL` site contributes an emission of 1 to every state — and is walked anyway.** Skipping
  it is not exact: the transition over two gaps is not the transition over their sum, because a
  state that leaves its class and returns redraws its templates, which one long step counts as kept.
  On the prototype, removing one flat site from a nine-site window moved the other sites' posteriors
  by up to 0.0005 (`S` = `R` = 30 kb). So every record with at least two rows is a site of every
  sample's walk, a sample's flat and absent records included (an absent sample's likelihood is flat
  too); the sample's posterior is kept only where it has a row. *Corrected 2026-11-01: this line
  said skipping was exact.*
- **Positions are unsigned.** A distance is computed from non-decreasing `POS` within a contig; a
  record out of order is a fault to report, not a negative gap.
- **Each record is written by exactly one window's core** (§2.2); writing one twice, or not at all,
  at a boundary is what §8's invariants check directly.

---

## 8. How we know it works

1. **Against the Python reference.** The corrected prototype's `hmm.py` gives the same posteriors to
   1e-9 on windows with multi-allelic records (`src/neighbour_evidence/model.rs`'s tests, done), and
   on an STR-bearing simulated window.
2. **On the founder cohort** (`str_stutter_benchmark.md` §12), scored against `truth.vcf` and
   `ancestry.tsv`: founders 4, 16, 64; F 0, 0.5, 0.8, 0.95; samples 10, 50, 200; depth 1, 3, 10, 30.
   It sets `N_min`, `F_min`, `K`, the window sizes, and whether frozen `S` and `R` cost accuracy. The
   autozygous posterior is scored against the runs in `ancestry.tsv`. **Two cases the founder model
   is built to make**: 10 samples from 64 founders, where most haplotypes have no template; and a
   mixed cohort, inbred lines with outbred samples (F 0) among them.
3. **Quality, not only error counts.** `GQ` against the observed error rate, in bins, separately for
   samples with reads and the no-read calls of §4.2; and the count of true-reference genotypes called
   non-reference, which the prototype tripled (1,382 to 4,227 at the nine-at-3× condition, report
   §5.2) — the owner sets what rise is acceptable (§11 Q5).
4. **On tomato1 held out**, as the report measured: the nine deep accessions cut to 3× and 1×, their
   full-depth calls as truth; and, with no downsampling, **the two plants sequenced both deep and
   shallow** (SRR7279542 with SRR7279488, SRR7279514 with SRR7279529), whose deep calls are truth for
   the shallow run. The prototype's numbers are a reference, not a floor, given §3.7.
5. **The circularity** (§3.3): at 10 and 50 samples, the pass with templates' `GP` recomputed
   without the target, against the pass as specified.
6. **Cost and back-pressure** at 57 and 500 samples on tomato1's regions, extrapolated to 3,000 by
   §6's arithmetic checked against the 500; how long the calling workers waited on the pass.
7. **Invariants:**
   - every record in is a record out, in order, unless §4.1 drops it;
   - output byte-identical at 1, 4 and 16 threads;
   - a record without rows or not PASS comes out byte-identical;
   - a long gap or a contig change leaves the autozygous probability at its long-run share (done,
     `a_long_gap_returns_the_autozygous_state_to_its_share`);
   - a sample with no informative record in reach keeps `./.` (§4.2);
   - accuracy at core edges no worse than at core centres.
8. **Refusal:** the stand-alone command refuses a VCF whose header says the pass ran (§9 Q4).

---

## 9. Decisions — **ruled by the owner, 2026-11-01**

All four as recommended.

1. **A sample whose reads said nothing: write the pass's genotype, or keep `./.`?** **Ruled: write
   it.** A call there is the aim the owner set ("where possible a call at a site with little or
   no evidence"), its `GQ` is honest about how sure it is, and its flat `PL` lets anyone filter it
   out. *The review's narrowing of this to calls with evidence nearby is §4.2 and §11 Q1.*
2. **Mark the refined genotypes per sample?** **Ruled: no new per-sample field.** The header says
   the pass ran and `GP`'s declaration says whose posterior it is (§4.1).
3. **`--neighbour-evidence auto|on|off`, with `auto` deciding by §5's two conditions?** **Ruled:
   yes.**
4. **A stand-alone command over a finished VCF** (`pvc refine-genotypes --vcf in.vcf.gz --parameters
   run.parameters.toml --output out.vcf.gz`)? **Ruled: yes, second, after the in-run stage.** It reads
   each sample's inbreeding coefficient from the run's parameters file, which the VCF names
   (`##parametersFile=`), and **refuses a VCF the pass has already refined**: refined `GP` as templates
   would count the cohort's evidence twice.

   *Built 2026-11-01 (`src/cli/refine_genotypes.rs`).* `--parameters` is required rather than read
   off `##pvcParameters=`, which records a path relative to wherever the run started; a file whose
   checksum differs from the one that line records is warned about, not refused, since the
   coefficients only start each sample's fit. The output keeps the input's header, with the
   `##neighbourEvidence=` line replaced, `GP`'s declaration changed and the command recorded in
   `##refineGenotypesCommandline=`. The pass always runs (it is the command's purpose); a cohort
   it cannot run on — one sample, a ploidy other than two — is an error naming why. On the
   20-sample simulated cohort its records are byte-identical to those the in-run pass wrote.

   *Removed 2026-11-03* (owner; read phase spec §5): refinement runs only inside a calling run,
   because read evidence never reaches a VCF, so a refinement of a file could not use it.

---

## 10. Deferred, with a recommended home

- **Template choice at thousands of samples by a positional Burrows–Wheeler index** (07 §2.6), and
  re-choosing templates more often than once a window, as the tools do every 0.1 cM — here, once §8
  measures the sampled choice at 500 samples.
- **Phased templates and sweeps over the cohort**, the outbred case — not this pass; 07 §2.
- **Ploidy above two** — a state per multiset of templates; here, when a polyploid cohort is a target.
- **Refining with the pass's own output** — a second sweep — here, if §8 shows the first leaves
  accuracy on the table; it must then be one controlled sweep, not a rerun (§9 Q4).

---

## 11. Open questions — **ruled by the owner, 2026-11-01**

Each with the recommendation it was ruled on. Q1–Q3 and Q5–Q8 as recommended; Q4: measure.

1. **Calls at no-read sites only with evidence nearby (§4.2)?** **Ruled: yes.** Without it the pass
   writes `0/0` at `GQ` about 17 where a sample has no informative record within reach — a failed
   library, a deletion in that line, a window at zero reads — which is the cohort's frequency, not
   evidence. It keeps every call the prototype recovered from neighbouring evidence.
2. **After the pass: drop alternatives no sample carries, and records with none left (§4.1)?**
   **Ruled: yes**, the rule calling already applies. The alternative is a record with `AC` 0 and a
   variant `QUAL`, which no reader can interpret.
3. **Mixed cohorts — ruled: run-level now; the mixed-cohort measurement of §8.2 decides whether to
   switch per sample.** A run-level switch by the median inbreeding runs the pass on a few outbred
   members of an inbred panel, where heterozygous templates blur, and turns it off for inbred lines
   in an outbred-majority cohort. *Recommend: run-level now, and the mixed-cohort measurement of §8.2
   decides whether to switch per sample* (each sample re-called only if its own fitted inbreeding is
   at least `F_min`, every sample still a template).
4. **Templates inherit calling's bias against heterozygotes.** A template's shares come from calling's
   `GP`, made under a genome-wide `F` near 0.95, which called 43% of true heterozygotes with four or
   more reads homozygous at 3× (report §5.2); a sample copying it inherits the error, and a template
   with no reads has shares near the frequency and can look like a good match. **Ruled: measure
   first** (§8.2 scores templates' own error); if it bites, compute template shares with the report's
   two-state autozygosity model before the walk, and weight agreement by the template's own
   informative records.
5. **What rise in false non-reference calls is acceptable?** With the nine accessions at 3×, the
   prototype fixed about 25,500 wrong genotypes at true non-reference sites and made about 2,800 true
   reference genotypes non-reference (report §5.2). **Ruled: judge by the net, measured per depth, and refuse a
   default that adds more false non-reference calls than it removes errors in any depth class.**
   **Restated by the owner, 2026-11-02:** a `./.` is not an error, and it is not a good outcome
   either. The aim is as little missing data as possible without filling it with wrong genotypes,
   so the pass is judged per depth by two numbers side by side: the share of `./.` it fills, and
   how many of the genotypes it fills are wrong — with, where samples have reads, the errors it
   removes against those it adds.
6. **No template carries the haplotype** — an introgression, a private lineage, 10 samples from 64
   founders. The mismatch rate is the only escape, so the pass imputes the nearest template's allele
   with confidence (07 §2.7 cites 18% heterozygous error in under-represented ancient genomes).
   **Ruled: add one background template whose shares are the cohort's allele frequencies** (§3.3),
   so a stretch no sample matches falls back to calling's own prior, **and confirm on the 64-founder
   case of §8.2.**
7. **Neighbouring records share reads.** At 26 records a kilobase a 150-base read covers about four;
   the walk treats their likelihoods as independent, so one mismapped read counts several times and
   `GQ` can come out too high at low depth. **Ruled: let §8.3's calibration decide**; if `GQ` is
   overconfident, thin the sites the walk uses for its state posterior to one per read length.
8. **The whole cohort at 1×.** The prototype was worse than the autozygous state alone at sites with
   one to three reads (22,408 against 20,111 wrong, report §5.2). **Ruled: measure on §8.2's 1×
   cases first**; Q4's template correction and Q6's background template are the likely cures, and a
   depth-dependent fallback is the last resort.
