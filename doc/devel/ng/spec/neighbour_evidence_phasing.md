# ng — phased templates for the neighbour-evidence pass

*Spec, 2026-11-02. No Rust yet — this settles the design. It extends
[`neighbour_evidence_pass.md`](neighbour_evidence_pass.md) (the pass) and rests on the Python
prototype [`phaser.py`](../research/neighbour_evidence/phaser.py) and the measurements in
[`phasing_for_the_pass_2026-11-02.md`](../research/neighbour_evidence/phasing_for_the_pass_2026-11-02.md).
The owner asked for it the same day, for partly outbred material — *S. pimpinellifolium* and
*cerasiforme* accessions among thousands of tomato accessions, many at medium to high coverage and
many at 5×. That narrows the ruling of 2026-09-13
([`07_decision`](../research/phasing_imputation/07_decision_separate_tool_and_route.md)) that phasing
stays outside pvc: phasing enters the pass as a means to better genotypes, not as an output.*

---

## 1. What this is

The pass re-calls each sample's genotypes by asking, along the chromosome, which other samples it
is copying. A template is now another sample's *genotype*. In an inbred sample a genotype is a
haplotype counted twice, so the copy is sharp. In a heterozygous sample it is not: a copy following
a template that is heterozygous at a site is told only "either allele, half and half". Partly outbred
material is heterozygous over much of its genome, and there the pass loses most of what it could
know.

**This adds phasing inside the pass:** each sample carries two haplotypes, its templates become the
haplotypes of the samples it resembles, and the haplotypes are estimated by rounds of the pass's own
walk. Nothing about the walk changes; what a template carries does.

**What the prototype measured** (Python reference model, 50 samples, founder simulator; wrong
genotypes where samples have reads):

| cohort | calling | pass as now | **prototype phaser** | Beagle 5.5 | perfect phasing |
|---|---|---|---|---|---|
| 10 at 30× + 40 at 3×, inbreeding 0.5, 16 founders | 2,187 | 1,106 | **699** | 787 | 392 |
| same, inbreeding 0.8 | 1,053 | 409 | **305** | 310 | 251 |
| same, inbreeding 0.8, 64 founders | 1,080 | 729 | **707** | 715 | 617 |
| all 3×, inbreeding 0.5, 16 founders | 2,595 | 1,339 | **906** | 1,154 | 472 |
| all 3×, inbreeding 0.8, 16 founders | 1,245 | 515 | **433** | 458 | 320 |
| all 3×, outbred, 16 founders | 4,934 | 3,971 | **3,352** | 3,755 | 591 |
| all 3×, inbreeding 0.95, 16 founders | 281 | 121 | 113 | 115 | 109 |
| all 3×, outbred, 64 founders | 5,308 | 4,934 | **5,069** | 5,263 | 3,003 |

And the `./.` it fills, first cohort: 1,589 of 2,116 at 0.8% wrong; the pass as now 1,260 at 1.5%,
Beagle's haplotypes 1,636 at 1.9%. The model's lengths are fixed at the simulator's and the pass's
rules for leaving a sample `./.` are not applied, so the numbers compare columns; they are not the
pass's own. The 64-founder cohorts give each haplotype one or two carriers among 50 samples; the
owner's point (2026-11-02) is that what matters is carriers per haplotype, which thousands of
accessions supply — simulations of 200 samples from 64 founders, 40 of them at 30×, are running
(§8).

### 1.1 Goals

- Give the pass haplotype templates wherever samples are heterozygous, estimated inside the pass
  from what it already reads (`PL`, `GP`, allele frequencies, each sample's inbreeding).
- Keep each haplotype allele's uncertainty: a copy is a probability per allele, not a hard allele.
  This is what makes the fills more reliable than Beagle's in the table above.
- Change nothing where it cannot help: an inbred cohort's result must stay as it is (row 7).
- Deterministic output, independent of thread count, as the pass is now.

### 1.2 Non-goals

- **Writing phase.** `GT` stays unphased and no `PS` is written: the haplotypes serve the genotype
  call. Writing them is *Deferred* (§9).
- **Phase across windows.** Each window phases its own records; nothing links one window's
  haplotypes to the next's, because nothing downstream needs it while phase is not written.
- **A reference panel.** The cohort phases itself, as the pass's templates come from the cohort.
- **Read-backed phasing** — reads spanning two heterozygous sites. Out, as for the pass (ruled
  2026-11-01).

### 1.3 Vocabulary

- **Haplotype** of a sample in a window: one of its two chromosome copies, as, at each record, a
  probability for each allele. Called *soft*: an allele can be uncertain.
- **Haplotype template**: one haplotype of another sample, used as a template exactly as a genotype
  template is used now — its allele shares at a record are the haplotype's probabilities.
- **Round**: one walk of every sample of the window, each against the haplotypes the previous round
  left, followed by every sample's haplotypes re-estimated.
- **Path**: one state per record, drawn at random from the walk — which template each copy follows
  at each record, or the one template both follow where the sample is autozygous.

---

## 2. The method

### 2.1 In words

1. **Start from what the pass does now.** Each sample's two haplotypes start as its allele shares
   from `GP`, the same for both copies. Round 1 is therefore the pass as now, with each template
   counted twice.
2. **Each round, for each sample:**
   - its templates are the two haplotypes of each of its `K` template samples (chosen once per
     window, as now, §2.3), plus the background template — `2K + 1` templates;
   - the walk gives its genotype posterior at every record, as now;
   - **one path is drawn** through the walk's states, backwards from the last record, each state
     drawn given the one after it and the forward probabilities;
   - **its two haplotypes are re-estimated along that path:** at each record, the copies follow
     templates `k₁` and `k₂` (or one template, autozygous), and the probability that copy one
     carries `a` and copy two `b` is the followed templates' probabilities times the sample's own
     likelihood of the genotype `{a, b}`, normalised:

     `P(a, b) ∝ q̃_{k₁}(a) · q̃_{k₂}(b) · L({a, b})`, and in an autozygous state
     `P(a, a) ∝ q̃_k(a) · L({a, a})`, `P(a, b) = 0` for `a ≠ b`,

     with `q̃` the templates' probabilities blurred by the mismatch rate `ε` as the walk blurs them
     (pass §3.4). Copy one's haplotype at the record is `P(a, ·)`, copy two's `P(·, b)`.
3. **All samples are updated together** at the end of a round: every walk of round `r` reads the
   haplotypes of round `r − 1`.
4. **After `R` rounds the genotype is the posterior averaged over rounds 2 to `R`.** Everything the
   pass writes and reports is taken from the averaged quantities (§4).

### 2.2 Decisions, with what they beat

- **A drawn path, not the posterior's marginals — decided.** The two copies are symmetric in the
  walk: averaged over states, copy one's allele probability equals copy two's at every record, and
  the phase is gone. One path breaks the symmetry; it is what SHAPEIT and Beagle sample too. The
  path is drawn, not the single most likely (Viterbi): a drawn path lets successive rounds explore,
  and the prototype was built and measured that way; Viterbi was not measured.
- **Soft haplotypes along the path, not hard alleles — decided, measured.** Beagle hands back hard
  haplotypes and its fills were wrong at 1.9% against the prototype's 0.8% in the same cohort; a
  copy that stays uncertain where neither its template nor the reads settle it does not pass a guess
  on as a certainty.
- **Average the posterior over rounds 2 to `R`, not take the last — decided, measured.** Samples
  phase each other, so an error can reinforce itself round to round: on the first cohort the last
  round's fills drifted from 1.4% wrong (round 3) to 2.8% (round 8), while the average stayed at
  0.8–1.1%. Round 1 is left out because it is the pass without phasing.
- **`R` = 5: a starting value**, from one cohort: the average over rounds 2–5 was 699 wrong, 2–7
  696, 2–8 712. Not measured elsewhere.
- **Round 1 from `GP` shares, both copies the same — decided.** The prototype's alternatives were
  the subtraction (one copy follows relatives, the other is what the genotype leaves) and Beagle;
  both were measured as *templates* and lost, neither was tried as a *start*. Starting from the
  pass as now makes round 1 a check: it must reproduce the pass's current output.

### 2.3 Templates

**The `K` template samples are chosen once per window as now** (pass §3.3: least disagreement where
both have reads, a pessimistic score), from `GP`, and each contributes both its haplotypes:
`2K + 1` templates with the background, 33 at `K` = 16. The prototype ran exactly this.

**Choosing 16 haplotypes instead of 16 samples is being measured** (§10 Q2, ruled: adopt it if within
a tenth): taking the 16 best of all `2(N − 1)` haplotypes keeps the walk's cost as it is now, and
needs a haplotype-to-genotype distance the code does not have.

### 2.4 What does not change

The walk itself (pass §3.4–§3.6: emission, transitions, posterior, the zero-emission guard), the
parameter fit (pass §3.5 — fitted in round 1, against genotype templates, and reused in later rounds
unchanged), windows and their sizes (§5 below changes the budget), the rewriting of the line, the
no-read rules (pass §4.2), and the fill check (pass §4.3).

---

## 3. Where it sits in the code

Inside the window solve (`refine_window`, `src/neighbour_evidence/refine.rs`). Today it solves each
sample once (`refine_sample`). With phasing it runs `R` rounds; each round solves every sample
against the window's current haplotypes, in parallel over samples as now
(`fits.par_iter_mut()`, `refine.rs:189`), and ends with a barrier where the round's haplotypes
replace the previous round's.

| what | existing code | how it is reused |
|---|---|---|
| the walk | `model::walk` (`model.rs:202`) | unchanged; called with `2K + 1` templates whose shares are haplotype probabilities |
| the forward probabilities a path is drawn from | `forward_autozygous`, `forward_pairs` inside `walk` (`model.rs:226–274`), dropped at its end | **must be returned** (or the draw done inside `walk`) — the one change to `model.rs` |
| the transitions a backward draw needs | the forward step's arithmetic (`model.rs:236–270`, `persistence`) | the draw uses the same transition, read backwards: given the next state, each current state's weight is its forward probability times the transition into the next |
| template choice | `templates::choose_templates` | unchanged: samples chosen from `GP`, each expanded to two haplotypes |
| a template's shares at a record | `SiteRecord::template_shares` (`records.rs:184`), from `GP` | round 1 only; later rounds read the window's haplotype store (§5) |
| the background template | `SiteRecord::background_shares` | unchanged, last |
| the genotype pairs of a record | `model::GenotypePairs` | the re-estimation's `L({a, b})` reads the sample's likelihood at the genotype index of `{a, b}` |
| the per-sample fit | `fit::SampleFit::fit_window` | round 1 only |
| the fill check, background-copy limit, no-read rules | `refine_sample` (`refine.rs`) | applied once, to the averaged prediction, posterior and background-copy probability (§4) |

The Python prototype is the oracle: `phaser.py`'s `draw_path` (the backward draw and the
re-estimation) and its round loop.

---

## 4. What the pass writes, with phasing

The averaged quantities replace the single walk's, everywhere the pass uses them:

- **the genotype posterior** → `GT`, `GQ`, `GP` (pass §4.1);
- **the neighbours-only prediction** → the fill check's comparison (pass §4.3);
- **the background-copy probability** → the no-read limit of 0.1 (pass §4.2). With haplotype
  templates "copying the background" keeps its meaning: no sampled haplotype matches.

Each is the mean of rounds 2 to `R`. A no-read call's `GQ` floor, the evidence-nearby rule and the
background limit are applied once, to the means. The header's `##neighbourEvidence=` line gains
`rounds=R`.

---

## 5. Cost, memory, determinism

**Time.** A walk's cost grows with the square of its templates: 33² against 17² states, 3.8 times.
Drawing a path is one more backward sweep, about the cost of the backward walk; the re-estimation
is `A²` a record. So a round is about 4.5 of today's walks and five rounds about 23. The pass's own
estimate at 3,000 samples is 21 core-hours (pass §6); with phasing about 480 — on the run's 30
threads, 16 hours. The owner (2026-11-02): parallelising over segments of the genome answers it.
Today windows are solved one after another inside the writer and samples in parallel within a
window (pass §2.3), which keeps every thread busy once a window has as many samples as threads;
windows solved in parallel are the next step if a cohort is smaller than that (§9). **The cost
estimate is arithmetic, not a measurement** — the first Rust build measures it (§8).

**Memory.**
- *The haplotype store*: two haplotypes a sample a record, `A` probabilities each, as `f32` — 16
  bytes a biallelic sample-record, against the 40 a window holds today (pass §6). Two copies of it
  live at once (the round being read and the round being written): 32 bytes, so a window holds
  about 72 bytes a biallelic sample-record. **The window budget (`budget_sample_records`,
  `windows.rs`) must fall to about 55% of its value to keep the same memory**, about 66 million
  sample-records; at 3,000 samples the window reaches its 500-kb floor sooner.
- *A walk's forward vectors*, kept for the backward walk and now for the draw: `T² + T` numbers a
  record — 1,122 at 33 templates. As `f64` over a 40,000-record window, **360 MB a sample, one per
  thread: 10.8 GB at 30 threads.** Today's 17 templates need 98 MB. Storing them as `f32` halves it;
  keeping only every `√n`-th record's vector and recomputing between them (checkpointing) brings it
  to megabytes at about one more forward sweep. **Which is a decision for the arch doc; the spec's
  requirement is that a thread's walk fits in a few hundred MB at 30 threads.**

**Determinism.** The draw is random. Each sample's path in each round is drawn from a generator
seeded by the run's fixed seed, the sample's index, the round, and the window's first core record —
nothing that depends on threads or on solving order. Rounds are synchronous (§2.1 step 3), so the
result is the same at 1 and 30 threads. **The seed is a constant in the code, not an option**; the
output is a function of the input.

---

## 6. When it runs

- **Only when the pass runs** (pass §5).
- **Always, with no option of its own — ruled** (owner, 2026-11-02, §10 Q1): "if we commit to the
  phasing model we'll use that model", and the command line stays as it is. An inbred cohort pays
  the cost for a small gain (row 7: 121 against 113 wrong); performance is to be worked on once the
  model runs, not by switching it off.
- **Diploid only**, as the pass.

---

## 7. Traps

- **The walk must not see a template twice by accident.** An inbred template's two haplotypes are
  near-identical; the walk treats them as two templates, which is correct (the prototype did) but
  doubles their prior weight against an outbred template's two distinct haplotypes. That is the
  prototype's behaviour and what was measured; do not "fix" it by merging without measuring (§9).
- **The draw reads the forward probabilities with the transition reversed.** Given the state at
  `i + 1`, the weight of state `s` at `i` is `forward_i(s) · P(s → state_{i+1})`. For a pair state
  `(k₁, k₂)` that transition has three terms (both copies keep, one keeps, neither), plus the
  class-change term from the autozygous states; `phaser.py`'s `draw_path` is the reference.
  A draw that uses the forward step's *margins* instead of the reversed transition draws a path
  that is not a path of the model.
- **The re-estimation needs the sample's own likelihood, and the posterior must not count it
  twice.** The re-estimated haplotypes contain the sample's reads; when *another* sample copies
  them, that is evidence about the other sample. When the sample itself is re-called its own
  haplotypes are never its templates (a sample is never its own template, pass §3.3) — keep it so.
- **But the circularity is larger than the pass's.** Sample `s` phases with `t`'s haplotypes, and
  `t` with `s`'s, so `s`'s reads return to it through `t` in the next round. The averaging over
  rounds keeps it from compounding in the measurements; the prototype did not leave the target out,
  and the Beagle comparison did (one run per target). Measured as part of §8, not removed.
- **Round 1 must reproduce today's pass exactly**, record for record, with each template doubled —
  the walk with a template counted twice is not the same walk (prior weights differ), so "exactly"
  means *the prototype's round 1*, not today's Rust output. Check against `phaser.py`.
- **An absent template sample** (no column at a record) has shares `1/A` for both haplotypes, as
  now; a flat one (no reads) gets haplotypes from its own path like any other sample.
- **Margins are phased in every window that reads them**, from scratch: the left margin of window
  `w + 1` was the core of window `w`. Re-using window `w`'s haplotypes there is tempting and wrong
  for determinism unless done identically at every thread count; leave it (§9).

---

## 8. How we know it works

1. **The prototype's numbers, reproduced.** On the founder-grid cells above, the Rust pass with
   phasing gives wrong-genotype counts within a tenth of the prototype's averaged column, the
   prototype run with the Rust pass's own rules and fitted parameters (to be added to `phaser.py`).
2. **Unit tests:** the backward draw on a case where one template equals the sample (the path must
   follow it everywhere, with any seed); the re-estimation on hand-computed values (a homozygous
   template followed by copy one and a heterozygous likelihood gives copy two the other allele —
   the owner's subtraction, as a special case); round 1 against the prototype's round 1 to 1e-9.
3. **Invariance:** output byte-identical at 1, 4 and 30 threads.
4. **No harm on inbred cohorts:** inbreeding 0.95 cells change by no more than the prototype's (121
   → 113).
5. **The ratio the owner pointed to:** 200 samples from 64 founders, 40 at 30×, inbreeding 0.5 and
   0.8 (running) — the case between the 50-sample cohorts and the owner's thousands.
6. **Cost**, measured: time and peak memory against the pass without phasing, at 50, 200 and 500
   samples; extrapolated to 3,000.
7. **Tomato, held out:** the tomato hiding test (`tomato_fill_test.sh`) on kimura, with and without
   phasing, on a cohort including *pimpinellifolium* and *cerasiforme* accessions.
8. **The fill check stays calibrated** with phasing: its estimated share against the true fill error
   on the cells of point 1, as pass §4.3 measured.

---

## 9. Deferred, with a recommended home

- **Writing the phase** (`GT` with `|`, `PS`): the haplotypes exist after the pass; writing them
  needs phase continuity across windows (§1.2). Home: a follow-up spec once phasing accuracy is
  measured on tomato against the deep accessions.
- **Windows solved in parallel**, for cohorts smaller than the thread count: the stage would hand
  windows to workers and reorder their lines. Home: pass §2.1's "a thread of its own" change.
- **Merging near-identical haplotype templates** (an inbred template's two copies) to cut the state
  count: a cost optimisation that changes the prior weights. Home: after §8.6 measures the cost.
- **Viterbi instead of a drawn path**, and **starting from the subtraction**: unmeasured
  alternatives (§2.2). Home: the prototype, if §8 falls short.

---

## 10. Questions — ruled by the owner, 2026-11-02

- **Q1. When phasing runs** (§6). *Ruled: always, whenever the pass runs, with no option.* The
  leaning offered — on when a tenth of the samples have inbreeding below 0.8 — was declined to keep
  the command line simple; performance is to be optimised instead.
- **Q2. 16 sample templates (32 haplotypes) or 16 haplotype templates** (§2.3). *Ruled: measure the
  second in the prototype and adopt it if its wrong-genotype counts are within a tenth of the
  first's.* It costs as now; the first costs 3.8 times a walk.
- **Q3. `R`** (§2.2). *Ruled: keep 5; measure 3 and 8 on the 200-sample cells.*
