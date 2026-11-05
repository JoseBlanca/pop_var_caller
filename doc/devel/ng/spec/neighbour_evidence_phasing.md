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

| cohort | calling | pass as now | **phaser, 16 haplotype templates** (adopted) | phaser, 32 | Beagle 5.5 | perfect phasing |
|---|---|---|---|---|---|---|
| 10 at 30× + 40 at 3×, inbreeding 0.5, 16 founders | 2,187 | 1,106 | **689** | 699 | 787 | 392 |
| same, inbreeding 0.8 | 1,053 | 409 | **322** | 305 | 310 | 251 |
| same, inbreeding 0.8, 64 founders | 1,080 | 729 | **695** | 707 | 715 | 617 |
| all 3×, inbreeding 0.5, 16 founders | 2,595 | 1,339 | **865** | 906 | 1,154 | 472 |
| all 3×, inbreeding 0.8, 16 founders | 1,245 | 515 | **369** | 433 | 458 | 320 |
| all 3×, outbred, 16 founders | 4,934 | 3,971 | **3,126** | 3,352 | 3,755 | 591 |
| all 3×, inbreeding 0.95, 16 founders | 281 | 121 | 115 | 113 | 115 | 109 |
| all 3×, outbred, 64 founders | 5,308 | 4,934 | 4,933 | 5,069 | 5,263 | 3,003 |

The `./.` it fills, first cohort: 1,485 of 2,116 at 0.8% wrong; with 32 templates 1,589 at 0.8%; the
pass as now 1,260 at 1.5%; Beagle's haplotypes 1,636 at 1.9%. In the outbred cohorts the 16
templates fill fewer and far more reliably: 764 at 1.7% wrong against 953 at 4.8% (16 founders),
579 at 5.4% against 733 at 8.7% (64). The model's lengths are fixed at the simulator's and the pass's
rules for leaving a sample `./.` are not applied, so the numbers compare columns; they are not the
pass's own. **They are biallelic SNP records only** (`phasing_ceiling.py` reads nothing else); the
pass phases every record, indels and multi-allelic STR records included, and those are unmeasured
(§8.9). The draw is random, so a cell moves with the seed: the first cohort gave 689, 724 and 687
wrong at three seeds, a spread of about 5% (review I7). The adopted column is
`uv run --with numpy python phaser.py --rounds 5 --haplotypes 16 --one-per-sample <cell>` on
`tmp/founder_grid/<cell>` of the worktree. The 64-founder cohorts give each haplotype one or two carriers among 50 samples; the
owner's point (2026-11-02) is that what matters is carriers per haplotype, which thousands of
accessions supply — simulations of 200 samples from 64 founders, 40 of them at 30×, are running
(§8).

### 1.1 Goals

- Give the pass haplotype templates wherever samples are heterozygous, estimated inside the pass
  from what it already reads (`PL`, `GP`, allele frequencies, each sample's inbreeding).
- Keep each haplotype allele's uncertainty: a copy is a probability per allele, not a hard allele.
  The fills' greater reliability than Beagle's in the table above is consistent with this; no arm
  isolates it (soft against hard haplotypes in the same phaser).
- Change nothing where it cannot help: an inbred cohort's result must stay as it is (row 7).
- Deterministic output, independent of thread count, as the pass is now.

### 1.2 Non-goals

- **Writing phase — in this spec.** The haplotypes serve the genotype call here; writing their phase,
  in confident blocks joined across windows, is
  [`neighbour_evidence_phase_output.md`](neighbour_evidence_phase_output.md) (owner, 2026-11-02).
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

1. **Round 1 is the pass as now**, exactly: each sample's `K` template samples (chosen as now, pass
   §3.3) contribute their `GP` allele shares. Each sample's two haplotypes then start from its
   path, as in every round.
2. **Each round, for each sample:**
   - from round 2 its templates are `K` *haplotypes* of other samples — the `K` whose alleles its
     genotype lacks least, at most one of each sample (§2.3) — plus the background template, so
     the walk has as many templates as the pass has now;
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
- **Soft haplotypes along the path, not hard alleles — decided.** A copy that stays uncertain where
  neither its template nor the reads settle it does not pass a guess on as a certainty. Consistent
  with Beagle's hard haplotypes filling at 1.9% wrong against the prototype's 0.8% in the same
  cohort, but Beagle differs in more than hardness. *Measured 2026-11-05:* reading hard alleles
  everywhere from round 2 cost the outbred 1× samples of a 40-sample cohort 12 homozygous-reference
  genotypes in 100 (72 → 60 right); the stored haplotypes stay soft.
- **From round 2 the states are weighed by definite alleles — decided 2026-11-05, measured.** A
  template's shares do two jobs in a walk: in the emissions they decide *whom the sample is
  copying*, and in the prediction and the drawn haplotypes they say *what it then carries*. Soft
  shares are right for the second and wrong for the first: an unsure template — a 1× sample at 50:50
  where it has no read — is never plainly wrong, so it explains the sample's reads as well as a
  template that truly carries them, and the walk spreads its weight onto it. So from round 2 each
  emission reads one allele drawn from each template's shares, the background's included, while the
  prediction and the drawn haplotypes read the shares themselves (§2.3). On the simulated F2 (two
  30× parents, 20 offspring at 1×) heterozygous offspring genotypes right went from 51 to 70 in 100;
  in two 40-sample natural cohorts with samples down to 1×, the 1× samples gained 2 to 4 genotypes
  in 100 more and the samples at 5× and above moved by at most 1 in 100. Drawing the background's
  allele too matters: left soft, it was the one template never plainly wrong, the walk copied it
  more, and the no-read calls it withholds (pass §4.2) rose from 8 to 21 in 100 of the outbred 1×
  homozygous-reference genotypes (`doc/devel/reports/implementations/refine_sources_2026-11-05.md`).
- **Average the posterior over rounds 2 to `R`, not take the last — decided, measured.** Samples
  phase each other, so an error can reinforce itself round to round: on the first cohort the last
  round's fills drifted from 1.4% wrong (round 3) to 2.8% (round 8), while the average stayed at
  0.8–1.1% (32 templates; with 16, 2.1% for round 5 alone — 35 of 1,645 — against 0.8% averaged).
  Round 1 is left out because it is the pass without phasing.
- **`R` = 5: a starting value**, from one cohort with 32 templates: the average over rounds 2–5 was
  699 wrong, 2–7 696, 2–8 712. Rounds 3 and 8 are to be measured on the 200-sample cells (§10 Q3).
- **Round 1 is the pass as now — decided.** The prototype's alternatives were the subtraction (one
  copy follows relatives, the other is what the genotype leaves) and Beagle; both were measured as
  *templates* and lost, neither was tried as a *start*. Starting from the pass as now makes round 1
  a check: without phasing (`R` = 1) its posteriors must be the current pass's, record for record.
  With phasing, round 1 differs from today's pass at the records hidden for the fill check (§4.1),
  where it reads a flat likelihood — and only there.

### 2.3 Templates

**Round 1: `K` = 16 template samples, as the pass chooses them now** (pass §3.3), each contributing
its `GP` allele shares.

**From round 2: 16 haplotypes, re-chosen every round, at most one of each sample — decided, measured
(§10 Q2).** A haplotype's score against the sample is a disagreement rate: at each record, the
probability that the haplotype carries an allele the sample's genotype lacks — for a biallelic
record `x · P(0/0) + (1 − x) · P(1/1)`, `x` the haplotype's alternative probability and the genotype
the sample's posterior from the previous round; in general `Σ_a h(a) · P(genotype without a)`. It
is summed over the sample's informative records (at most 500, spread evenly, as now) where the
haplotype's sample is informative too, and scored pessimistically as the pass scores samples (the
Beta mean plus two standard deviations). Each candidate sample offers its better-scoring haplotype;
the 16 best samples' offers are the templates.

**What a round's walk reads of its templates.** The emissions read *definite alleles*: at each site,
one allele drawn from each template's shares — the 16 haplotypes' and the background's — from the
stream `StreamPurpose::DefiniteAlleles`, keyed by the window, the sample and the round (§5), one
uniform a template a site in template order. The genotype prediction and the re-estimation of the
sample's own haplotypes read the shares (§2.2). Round 1 reads the shares for both: its templates are
genotypes, not haplotypes, and it is the pass without phasing.

**A sample not walked in the window** — no informative record there, so no templates and no path
(`templates.rs:46`, `refine.rs:243`) — **carries its `GP` allele shares as both haplotypes in every
round**, the prototype's starting value, and is a candidate like any other. It scores worst (no
shared records), but in a cohort of 17 samples or fewer every other sample is chosen.

- *Why at most one of each sample:* without it the 16 slots fill with copies of the commonest
  haplotype — many inbred samples carry it twice — and the variety a heterozygous sample needs is
  lost: 871 wrong against 689 on the first cohort of §1.
- *Why not both haplotypes of 16 samples:* measured (the "32" column of §1), within a tenth of the
  adopted rule everywhere, but a walk with 33 templates costs 3.8 times one with 17.

### 2.4 What does not change

The walk itself (pass §3.4–§3.6: emission, transitions, posterior, the zero-emission guard), the
parameter fit (pass §3.5 — fitted in round 1, against genotype templates; later rounds read
`SampleFit::parameters` and must not call `fit_window` again, which would add the window's
log-likelihoods to its totals a second time, `fit.rs:122`), windows and their sizes (§5 below changes the budget), the rewriting of the line, the
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
| the walk | `model::walk` (`model.rs:202`) | unchanged; called with `K + 1` templates, from round 2 haplotype probabilities |
| the forward probabilities a path is drawn from | `forward_autozygous`, `forward_pairs` inside `walk` (`model.rs:226–274`), dropped at its end | **must be returned** (or the draw done inside `walk`) — the one change to `model.rs` |
| the transitions a backward draw needs | the forward step's arithmetic (`model.rs:236–270`, `persistence`) | the draw uses the same transition, read backwards: given the next state, each current state's weight is its forward probability times the transition into the next |
| template choice | `templates::choose_templates` | round 1, unchanged; from round 2 a new function choosing haplotypes (§2.3), reusing its record sampling and its pessimistic score (`templates::disagreement_bound`) |
| a template's shares at a record | `SiteRecord::template_shares` (`records.rs:184`), from `GP` | round 1 only; later rounds read the window's haplotype store (§5) |
| the background template | `SiteRecord::background_shares` | unchanged, last |
| the genotype pairs of a record | `model::GenotypePairs` | the re-estimation's `L({a, b})` reads the sample's likelihood at the genotype index of `{a, b}` |
| the per-sample fit | `fit::SampleFit::fit_window` | round 1 only |
| the fill check, background-copy limit, no-read rules | `refine_sample` (`refine.rs`) | applied once, to the averaged prediction, posterior and background-copy probability (§4) |

The Python prototype is the oracle: `phaser.py`'s `draw_path` (the backward draw and the
re-estimation) and its round loop.

### 3.1 The interfaces the change needs

- **`model::walk_and_draw(sites, templates, parameters, uniforms) -> (WalkResult, haplotypes)`**, a
  second function beside `walk`, sharing its forward code. The draw needs the forward vectors, the
  private `persistence` (`model.rs:183`) and the blurred shares in the private `SiteTerms`, so it
  belongs inside `model.rs`; `walk` itself stays as it is, so the fit's grid walks (`fit.rs:112–121`)
  do not carry vectors they never read. `uniforms` is the stream of §5's generator; the returned
  haplotypes are two vectors of `A` probabilities a walked record.
- **`refine_sample` splits in two.** *A round's solve*, per sample: its templates (round 1 as now,
  later the haplotype choice), `walk_and_draw`, and back come the new haplotypes at every walked
  record, the posterior at its compared records (for the next round's choice), and the posterior,
  neighbours-only prediction and background-copy probability at core records (to accumulate). *The
  final step*, per sample: the means over rounds 2 to `R`, then the no-read rules and the fill check
  as `refine_sample` applies them now, building `SampleSolution` as now.
- **The background's position** in a walk's shares is after the round's templates:
  `settled_genotype`'s slice (`refine.rs:325`) uses the round's template count, not `templates.len()`
  of round 1.
- **`R` is an internal constant, and `R` = 1 means "write round 1"**, so the existing `refine.rs`
  tests, run at `R` = 1, are §7's check that round 1 is today's pass.

---

## 4. What the pass writes, with phasing

The averaged quantities replace the single walk's, everywhere the pass uses them:

- **the genotype posterior** → `GT`, `GQ`, `GP` (pass §4.1);
- **the neighbours-only prediction** → the fill check (pass §4.3), **but only at hidden records**
  (§4.1): elsewhere, with phasing, the prediction is no longer free of the sample's own reads (§7);
- **the background-copy probability** → the no-read limit of 0.1 (pass §4.2). With haplotype
  templates "copying the background" keeps its meaning: no sampled haplotype matches.

Each is the mean of rounds 2 to `R`. A no-read call's `GQ` floor, the evidence-nearby rule and the
background limit are applied once, to the means; neither was applied in the prototype, so their
values are unmeasured with haplotype templates (§8.10). The header's `##neighbourEvidence=` line
gains `rounds=R` and `haplotype_templates=16`, so a reader can tell which model made the file.

### 4.1 The fill check, on hidden records

**Before round 1, up to 50 records of each sample in each window are hidden, and at most 5% of
those where its reads are informative** (owner, 2026-11-02, §10 Q4; the numbers set at the Rust
build): drawn at random among the window's core records where the sample's reads are informative, from
§5's generator. At a hidden record the sample is treated as having no reads in every round — its
walk and its path read a flat likelihood there — so its own walk never sees them, and its averaged
neighbours-only prediction there is close to what a real `./.` would be filled with: no reads at
the record, reads around it. *Close, not equal:* in round 1 the other samples copy its `GP` shares,
which calling made from those reads, and carry them into their round-1 haplotypes. Masking those
shares too was measured (code review 2026-11-02, 6,000 hidden records on the crate's test cohorts:
prediction errors 144 → 141, 34 → 43, 9 → 17) and showed no systematic effect, so the `GP` stays.

- **The fill check runs at hidden records only:** where the sample's own evidence settles its
  genotype and the §4.2 rules would write the prediction, the prediction's genotype against the
  settled one. Nothing else about pass §4.3 changes — the table, the log, the 30-site floor, the
  warnings.
- **The genotype written at a hidden record** is the averaged prediction times the sample's own
  likelihood there, renormalised — a posterior like any other, its reads entering once.
- **What it costs:** the sample's reads at a hidden record are missing from what other samples copy
  there. In tomato's windows 50 is about 1 in 1,000 of a sample's records. In the Rust pass on the
  2-Mb simulated cohorts (the 5% cap binding, about 25 a contig) wrong genotypes moved from 398 to
  401, 615 to 587 and 69 to 69 against 20 a window, and from 2,613 to 2,809 on the outbred cohort;
  in the prototype, hiding 200 of about 1,000 records a sample (a fifth) added about 30%.
- **Why 50 a window, capped, and not a fixed rate:** a rate small enough to be free over a genome (1
  in 1,000) gives a few-megabase run a handful of records a sample, and the check needs 30; the
  first build's 20 a window could never reach 30 in a run of one window. With 50 and the cap, a
  single-window run estimates for 12 to 28 of 50 samples on the inbred and partly inbred cohorts and
  for none on the outbred one, where few genotypes settle at 3×; a whole-genome run sums hundreds of
  windows.

Measured in the prototype (`phaser.py --hide 200`, hiding many for the measurement only), the
check's estimate against the fills' true error:

| cohort | true | check |
|---|---|---|
| 10 at 30× + 40 at 3×, inbreeding 0.5 | 1.1% | 0.6% |
| same, inbreeding 0.8 | 0.7% | 0.4% |
| same, 0.8, 64 founders | 2.5% | 1.5% |
| all 3×, inbreeding 0.5 | 1.4% | 0.8% |
| all 3×, inbreeding 0.8 | 0.6% | 0.6% |
| all 3×, inbreeding 0.95 | 0.5% | 0.4% |
| all 3×, outbred, 16 founders | 3.7% | 6.0% (133 records) |
| all 3×, outbred, 64 founders | 6.8% | 11.8% (170) |

Within a factor of two everywhere and in the right order: low where the partly inbred cohorts' fills
are good, high where the outbred ones are poor. The echoed check of §7 read 10.2% against 1.7% on
the outbred cohort. Its known biases: it reads low in partly inbred cohorts — a settled genotype at
3× is mostly a homozygote, the easier fill — and high in outbred ones; at 3× in outbred samples few
genotypes settle at all.

---

## 5. Cost, memory, determinism

**Time.** A walk has as many templates as now. Drawing a path is one more backward sweep, about a
third of a walk; the re-estimation is `A²` a record. So a round is about 1.3 of today's walks and
five rounds about 6.5. Choosing haplotypes each round adds `N · 2 · 500` comparisons a sample a
window: at 3,000 samples and the 1,800 windows of pass §6, about 5 core-hours over four rounds. At
that cohort size the windows are at their floor — a 0.5-Mb core in a 1.5-Mb window — so every record
is walked three times, not the two pass §6 counted: the pass costs about 33 core-hours, and with
phasing about 213 — on the run's 30 threads, about 7 hours. The owner (2026-11-02): parallelising
over segments of the genome answers what remains.
Today windows are solved one after another inside the writer and samples in parallel within a
window (pass §2.3), which keeps every thread busy once a window has as many samples as threads;
windows solved in parallel are the next step if a cohort is smaller than that (§9). **The cost
above is arithmetic; what the first Rust build measured** (code review, 2026-11-02, the crate's test
cohort on one thread): five rounds cost 3.8–5.6 times one round with the fit frozen, and 1.45 times
with round 1's parameter fit included (120 samples, 2,000 records: 23.0 s against 33.4 s) — the fit's
twelve grid walks dominate a first window. The haplotype choice first ran at about 30 ns a
comparison, against the 0.3 ns the arithmetic assumed, from two binary searches per comparison —
about 500 core-hours at 3,000 samples; the searches are now table lookups, not yet re-measured.

**Memory.** Three new holdings, and one the pass already had uncounted.
- *The haplotype store*: two haplotypes a sample a record, **`A − 1` probabilities each** (the last
  allele's is one minus the rest), as `f32` — 8 bytes a biallelic sample-record; two copies live at
  once (the round being read and the round being written), 16 bytes, against the 40 a window holds
  today (pass §6). **At absent rows the store holds `1/A`**, written when the store is filled, so a
  reader never has to ask; a sample's own walk overwrites only its walked records. *Built so*
  (owner, 2026-11-02, after a first build stored all `A` probabilities as `f64`, 64 bytes a biallelic
  sample-record for the two buffers): on the simulated cohorts the change moved one genotype of
  175,146 on the 200-sample cohort and none on three 50-sample ones. Each sample writes its own row
  of the round being written from its own thread, so no second copy of every sample's haplotypes is
  held at a round's end (the first build held one, about 80–110 bytes a sample-record; code review
  2026-11-02).
- *The round accumulators*: running sums, per sample and core record, of the posterior, the
  prediction and the background-copy probability — 7 numbers biallelic, as `f32` 28 bytes a core
  sample-record. The previous round's posterior is kept **only at each sample's at most 500
  compared records**, which is all the haplotype choice reads.
- *A walk*, per thread: the forward vectors (`T² + T` = 306 numbers a record at 17 templates) and
  the per-site terms the walk keeps for the backward sweep (`SiteTerms`, `model.rs:212–222`: 340
  numbers), as `f64` about 207 MB a sample over a 40,000-record window — **6.2 GB at 30 threads**.
  The draw reads both after the backward sweep, so they live longer, not larger. This was already so
  without phasing; pass §6 counted only the forward vectors.

At 3,000 samples the window is already at its floor (pass §6: about 39,000 records, 117 million
sample-records), so a smaller budget cannot shrink it: the window holds about 117 million × 56 bytes
≈ 6.6 GB, the accumulators about 1.1 GB over a 13,000-record core, and the walks 6.2 GB — about
14 GB, against about 4.7 + 6.2 GB without phasing. Below about 2,000 samples, lowering
`budget_sample_records` (`windows.rs`) to 70% of its value keeps a window's memory as it is now.

**Determinism.** The draw is random. Each sample's path in each round reads uniforms from a
SplitMix64 generator (the crate has no random-number dependency; it already uses hand-written
SplitMix64, e.g. `src/cli/estimate_contamination.rs:1005`) seeded by mixing a constant, **the
window's key — its contig's number and its first core site's position** (`refine::window_key`), the
sample's index and the round. The key depends on the window alone, so a window's draws do not
change with what else the run called; the first build keyed by the run-wide number of the window's
first site, which changed a contig's draws with the contigs called before it (code review,
2026-11-02). *The genotypes still can*: each sample's walk parameters are fitted on the run's first
windows and then frozen (pass §3.5), so calling a contig alone or after others fits them on
different windows. That predates phasing.
Each state is drawn with one uniform in [0, 1) against the cumulative weights in a
fixed order — autozygous states first, then pair states row by row — as `phaser.py` orders them.
The definite alleles a round's emissions read (§2.3) come from a stream of their own, keyed the same
way with its own purpose number, so drawing them changes no path.
Nothing depends on threads or on solving order, and rounds are synchronous (§2.1 step 3), so the
result is the same at 1 and 30 threads. **The seed is a constant in the code, not an option.**

---

## 6. When it runs

- **Only when the pass runs** (pass §5).
- **Always, with no option of its own — ruled** (owner, 2026-11-02, §10 Q1): "if we commit to the
  phasing model we'll use that model", and the command line stays as it is. An inbred cohort pays
  the cost for a small gain (row 7: 121 against 115 wrong); performance is to be worked on once the
  model runs, not by switching it off.
- **Diploid only**, as the pass.

---

## 7. Traps

- **The haplotype choice reads the previous round's posterior, not this round's**, and the
  previous round's haplotypes: every sample of round `r` chooses from round `r − 1`, so the choice
  is the same at any thread count and any order.
- **Ties in the haplotype score** break by sample index, then copy (0 before 1), as the prototype's
  sort does.
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
- **But the circularity is larger than the pass's, and it reaches the fill check.** A shallow sample
  `t` that follows `s`'s haplotype at a record re-estimates its own haplotype there as nearly `s`'s,
  which `s`'s reads shaped; next round `s` picks that haplotype, because it disagrees with `s`'s
  genotype least; and `s`'s neighbours-only prediction at that record echoes `s`'s own reads. For
  the genotype the averaging over rounds keeps it from compounding in the measurements (the
  prototype did not leave the target out; the Beagle comparison did). **For the fill check it
  breaks the premise** — a prediction without the sample's reads: on the first cohort a single
  phased round's true fill error stayed at 1.4% while the check's estimate fell to 0.65%; on the
  outbred 3× cohort the averaged check read 10.2% against a true 1.7% (review I1). §10 Q4.
- **Round 1 must reproduce today's pass exactly**, record for record: the same templates, the same
  walk. A difference there is a bug in the round loop, not in phasing.
- **An absent template sample** (no column at a record) reads `1/A` from the store (§5); a sample
  without reads at a record gets haplotypes there from its own path like any other record; a sample
  without a walk in the window carries its `GP` shares (§2.3).
- **Margins are phased in every window that reads them**, from scratch: the left margin of window
  `w + 1` was the core of window `w`. Re-using window `w`'s haplotypes there is tempting and wrong
  for determinism unless done identically at every thread count; leave it (§9).

---

## 8. How we know it works

1. **Exact agreement with the prototype — built, then retired.** `phaser_reference.py` uses §5's
   SplitMix64 stream and state order; on a drawn window of 20 samples and 500 records the Rust
   posteriors and hidden records agreed with it to 1e-9 through the first build and its code-review
   fixes (5f242a4a). The `f32` store (§5) moves the Rust results in their last digits, and the owner
   ruled the comparison no longer needed (2026-11-02): the test and its fixture are removed; the
   script stays as the method's record. Records with more than two alleles are covered by unit tests
   of the draw, the re-estimation and the haplotype choice. This, not a statistical tolerance, is what tells a bug from chance.
2. **The prototype's numbers, reproduced** on the founder-grid cells above, with the prototype run
   under the Rust pass's own rules and fitted parameters (to be added to `phaser.py`): wrong-genotype
   counts within a tenth, judged against the seed spread (about 5%, §1).
3. **Unit tests:**
   - the backward draw with `mismatch` = 0, one template equal to the sample and every other
     template and the background carrying the other allele: every drawn state is the autozygous
     state on that template or the pair state on it twice, at any seed;
   - the re-estimation on hand-computed values, the owner's subtraction as a special case: copy one
     on a homozygous-reference template, copy two on a template uniform over the alleles, and a
     likelihood allowing only the heterozygote give copy two the alternative with probability
     `(1 − ε/2)·½ / ((1 − ε/2)·½ + (ε/2)·½)` = 1 − ε/2 for copy one's reference — computed in the test;
   - round 1 against today's pass: the existing `refine.rs` tests at `R` = 1 (§3.1).
4. **Invariance:** output byte-identical at 1, 4 and 30 threads.
5. **No harm on inbred cohorts:** inbreeding 0.95 cells change by no more than the prototype's (121
   → 115).
6. **The ratio the owner pointed to:** 200 samples from 64 founders, 40 at 30×, inbreeding 0.5 and
   0.8 (running) — the case between the 50-sample cohorts and the owner's thousands.
7. **Cost**, measured: time and peak memory against the pass without phasing, at 50, 200 and 500
   samples; extrapolated to 3,000.
8. **Tomato, held out:** the tomato hiding test (`tomato_fill_test.sh`) on kimura, with and without
   phasing, on a cohort including *pimpinellifolium* and *cerasiforme* accessions.
9. **Records other than biallelic SNPs:** a multi-allelic case and an STR record through the haplotype
   choice's general formula and the re-estimation at `A` > 2 (unit), and the founder cells' indel
   and STR records scored on and off (the prototype never read them).
10. **The no-read rules with haplotype templates:** fills kept and their true error under the 0.1
    background limit and the `GQ` 10 floor, phasing on against off; and wrong calls per `GQ` bin on
    and off, since the sample's own reads now enter twice.
11. **Edge cohorts:** 2, 3 and 17 samples (every sample a template); absent columns; a sample with
    no informative record in a window (§2.3).
12. **Leave-the-target-out** on a few targets of one cell, against the run as built: the size of the
    circularity of §7 on the genotype.
13. **Rounds 3 and 8** against 5 on the 200-sample cells (§10 Q3).
14. **The fill check on hidden records** (§4.1): its estimate against the true fill error on the cells
    of point 2, the 200-sample cells included, at the production 50 a window — the prototype
    measured it with 200 a sample; and that hiding 50 a window moves wrong-genotype counts by less
    than the seed spread.

---

## 9. Deferred, with a recommended home

- **Writing the phase** (`GT` with `|`, `PS`): now
  [`neighbour_evidence_phase_output.md`](neighbour_evidence_phase_output.md).
- **Windows solved in parallel**, for cohorts smaller than the thread count: the stage would hand
  windows to workers and reorder their lines. Home: pass §2.1's "a thread of its own" change.
- **Viterbi instead of a drawn path**, and **starting from the subtraction**: unmeasured
  alternatives (§2.2). Home: the prototype, if §8 falls short.

---

## 10. Questions — ruled by the owner, 2026-11-02

- **Q1. When phasing runs** (§6). *Ruled: always, whenever the pass runs, with no option.* The
  leaning offered — on when a tenth of the samples have inbreeding below 0.8 — was declined to keep
  the command line simple; performance is to be optimised instead.
- **Q2. 16 sample templates (32 haplotypes) or 16 haplotype templates** (§2.3). *Ruled: measure the
  second and adopt it if within a tenth of the first.* **Resolved: 16 haplotype templates, at most
  one of each sample.** Within a tenth on all eight cohorts of §1 (on the first, 689 against 699
  wrong — a difference within the seed spread); fills fewer but more reliable. The 16 best haplotypes without the
  one-a-sample rule were 25% worse (871) and are not used.
- **Q3. `R`** (§2.2). *Ruled: keep 5; measure 3 and 8 on the 200-sample cells.*
- **Q4. What the fill check compares with phasing on** (§4.1). The review found the averaged
  prediction echoes the sample's own reads (§7). Offered: round 1's prediction, the averaged one, or
  removing the echo at a cost of a round. *The owner's proposal, adopted after the prototype
  measured it:* hide a few records of each sample from the start, so the echo never forms there,
  and check only those — the hiding test of Beagle users, built in. Up to 50 a window (at most 5%)
  rather than the proposed 1 in 1,000, so small runs get enough.
