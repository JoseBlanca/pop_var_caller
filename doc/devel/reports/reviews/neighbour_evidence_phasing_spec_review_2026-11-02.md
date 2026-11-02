# Review: phased templates for the neighbour-evidence pass (spec of 2026-11-02)

*Adversarial review of [`neighbour_evidence_phasing.md`](../../ng/spec/neighbour_evidence_phasing.md)
against the pass's code (`src/neighbour_evidence/`), the pass spec, and the prototype
(`doc/devel/ng/research/neighbour_evidence/phaser.py`, `phasing_ceiling.py`, `hmm.py`). The owner's
rulings (always on, 16 haplotype templates one per sample, `R` = 5, no phase written) are taken as
given. Probe scripts and their outputs are in `tmp/phasing_review/` of the `neighbour-pass` worktree.*

## Summary

**No blocker.** The method is described correctly and the prototype's arithmetic is sound:

- **The backward draw is a valid sample of the model's path.** I drew 30,000 paths with
  `phaser.draw_path` on small random cases (8 records, 3 templates plus the background) and averaged
  the genotype each drawn state implies; the average matched `hmm.haplotypes`' posterior to within
  0.0034 everywhere, against a Monte Carlo error of about 0.003 (`tmp/phasing_review/check_draw.py`).
  The forward step of `draw_path` (`phaser.py:75–84`) is the same arithmetic as
  `model.rs:236–260` and `hmm.py:185–197`, so a Rust port of the draw on `walk`'s forward vectors is
  correct. Gaps of zero (two records at one position) make the transition the identity, which the
  draw handles without special code.
- **Spec and prototype agree** on round 1, the haplotype choice and its tie order, the
  re-estimation formula (`phaser.py:120–130`) and the averaging over rounds 2 to `R`
  (`phaser.py:278–288`).
- **The headline number reproduces:** with `--rounds 5 --haplotypes 16 --one-per-sample` on the
  first cohort I get 689 wrong and 1,485 fills at 12 wrong (0.8%), as the spec says.

Nine findings need action before coding. The most important: phasing changes what the fill check
measures (I1), the background limit and `GQ` floor were never applied in the prototype (I2), the
memory and time figures in §5 are low by about a factor of two (I3, I4), and the spec does not
say what a sample that is never walked contributes as a template (I5).

---

## Important

### I1. Phasing feeds a sample's own reads into its "neighbours-only" prediction, so the fill check no longer measures what pass §4.3 says it does

**Section:** §4, §7 (circularity trap), §8.8.

**The mechanism.** The fill check compares, at records where the sample's reads settle its
genotype, the neighbours' prediction with that genotype (`refine.rs:321–335`). It is honest only if
the prediction contains none of the sample's own reads there. With phasing it does, and at the
record itself:

- in round `r − 1`, a shallow sample `t` that follows `s`'s haplotype at a record re-estimates its own
  haplotype there as ∝ `q̃_s(a) · L_t(…)`. With `t`'s likelihood nearly flat, that is `s`'s haplotype,
  which `s`'s reads shaped;
- in round `r`, `s` picks `t`'s haplotype as a template, because it disagrees with `s`'s genotype
  least (§2.3);
- `s`'s prediction at that record then echoes `s`'s own reads.

The spec's §7 discusses circularity only as a matter of accuracy and states that it is measured in
§8. **No item in §8 measures it**, and §4 does not say that the fill check is affected.

**Measured** (`tmp/phasing_review/leak_check.py`, prototype rounds, a fill check computed as
`refine.rs` computes it except for the evidence-nearby rule and the background limit). First
cohort (10 at 30× + 40 at 3×, inbreeding 0.5):

| round | true error of fills (records without reads) | fill check's estimate |
|---|---|---|
| 1 (pass as now) | 19 of 1,260 (1.5%) | 150 of 10,319 (1.5%) |
| 2 | 21 of 1,463 (1.4%) | 80 of 12,333 (0.65%) |
| 5 | 35 of 1,645 (2.1%) | 96 of 12,709 (0.76%) |
| mean of 2–5 | 12 of 1,485 (0.8%) | 104 of 12,035 (0.86%) |

In a single phased round the true fill error stays at 1.4% while the check's estimate halves, to
0.65%. That is the leak. Averaged over rounds 2–5 the two happen to agree on this cohort. On the
outbred 3× cohort they do not: the check gives 10.2% (81 of 798) against a true 1.7% (13 of 764).
My probe's absolute values differ from pass §4.3's because it applies neither the evidence-nearby
rule nor the background limit, so compare rows within a cohort, not with the pass spec. The
conclusion stands: **how the check relates to the true error depends on the cohort once phasing
is on.**

**Fix.** State the mechanism in §7 and choose a rule in §4. The cheap option is to compute the fill
check from round 1's prediction, which is the pass as now and calibrated as now. On the first
cohort it overstates the phased fills' error (1.5% against 0.8%), which errs on the safe side.
Alternatively keep the averaged prediction and make §8.8 a gate, with the outbred cells and the
200-sample cells, before the check's figure is reported with phasing on.

### I2. The background-copy limit and the `GQ` floor were never applied in the prototype, so their values are uncalibrated with haplotype templates

**Section:** §2.4, §4.

§2.4 lists the no-read rules as unchanged, and §4 applies them to averaged quantities. But the
prototype applies neither rule: `phaser.py` scores a fill at `GQ` ≥ 10 alone (`phaser.py:222`), with
no background limit and no evidence-nearby rule. The spec itself says so at line 46. The limit of
0.1 (`refine.rs:57`) was set on genotype templates, where each template stands for two
haplotypes. With 16 single haplotypes, a heterozygous sample is more likely to find no matching
haplotype and copy the background, so the same limit may withhold a different share of fills.
The averaged background-copy probability is a mixture over rounds — a legitimate probability, but
not one anyone has measured against the 0.1 threshold.

**Fix.** Add the two rules to `phaser.py`, which §8.1 already plans, and add an explicit §8 item:
the fills kept and their true error under the 0.1 limit and the `GQ` 10 floor, with phasing on
against off.

### I3. §5's memory figures are low: the walk holds twice what is stated, the budget cut cannot act at 3,000 samples, and the per-round accumulators are missing

**Section:** §5, *Memory*.

1. **The walk's own memory.** The spec counts the forward vectors, `T² + T` = 306 `f64` a record
   (98 MB a sample over 40,000 records, correct for the code). But `walk` also keeps every site's
   `SiteTerms` for the whole window (`model.rs:212–222`): the blurred shares, the `T` autozygous
   emissions and the `T²` pair emissions, 34 + 17 + 289 = 340 `f64`, 2,720 bytes a record — **as
   large as the forward vectors**. Together that is about 207 MB a sample at 40,000 records, so
   **about 6.2 GB at 30 threads, not 2.9 GB.** The re-estimation needs the blurred shares, so they
   must stay alive until the draw finishes.
2. **The budget cut does not act at 3,000 samples.** Pass §6 sizes windows at 3,000 samples to the
   floor already: 1,800 windows over 900 Mb are 0.5-Mb cores with 0.5-Mb margins, a 1.5-Mb window of
   about 39,000 records, which is 117 million sample-records against the budget of 120 million.
   Lowering the budget to 66 million (§5) cannot shrink the window further: the planner stops at the
   floor and reports the overrun (`windows.rs:124–125`). So at 3,000 samples the window holds
   117 million × 72 bytes ≈ 8.4 GB, against about 4.7 GB without phasing. The sentence "the window
   reaches its 500-kb floor sooner" is true below about 1,700 samples (66 million ÷ 39,000); at
   3,000 the cohort is past it.
3. **What the rounds accumulate is not counted.** Averaging needs, for every sample at every core
   record, the running sums of the posterior, the prediction and the background-copy probability:
   3 + 3 + 1 numbers biallelic, 56 bytes as `f64`. At 3,000 samples and a 13,000-record core that is
   about 2.2 GB. Today's `refine_window` already collects every sample's `SampleSolution` before
   assembling (`refine.rs:188–202`), about 56 bytes a core sample-record, and pass §6 does not count
   that either. The haplotype choice also needs the previous round's posterior, but only at each
   sample's at most 500 compared records — the spec should say so, or a coder will keep the whole
   window's posteriors (24 bytes a sample-record).

**Fix.** Restate §5's memory as: window ≈ (40 + 32) bytes a sample-record, plus about 56 bytes a
core sample-record in accumulators, plus about 210 MB a thread. Give the 3,000-sample total
(≈ 8.4 + 2.2 + 6.2 GB) and say that the budget cannot hold it there. Two cheap savings to name:
store each haplotype as `A − 1` numbers (biallelic: one `f32`, halving the store to 8 bytes,
16 double-buffered), and accumulate in `f32`.

### I4. The time estimate inherits pass §6's "every record walked twice"; at 3,000 samples it is three times

**Section:** §5, *Time*; the error is in pass §6 (`neighbour_evidence_pass.md:402`).

At the floor-sized windows that 1,800 windows imply, the window is 1.5 Mb and the core 0.5 Mb, so
each record is walked in three windows, not two. The pass's 21 core-hours becomes about 31, and
phasing's 140 becomes about 205, about 7 hours on 30 threads rather than under 5. The haplotype
choice's 5 core-hours is correct under the pass's own counting (3 × 10⁶ comparisons a sample a
window, × 3,000 samples × 1,800 windows × 4 rounds = 6.5 × 10¹³, at 4 × 10⁹ a second).

**Fix.** Correct both specs: three walks a record at the floor, two at the default sizes.

### I5. A sample with no informative record in the window is never walked, so it has no haplotypes — yet it is a candidate template

**Section:** §7 (absent and flat templates), §2.3.

`choose_templates` returns `None` for such a sample (`templates.rs:46–47`), and `refine_sample` then
returns `None` (`refine.rs:243`). There is no walk, so no path and no haplotypes. The prototype
never meets this case: `pc.choose` always returns templates, and every sample is walked
(`phaser.py:183–204`). §7's sentence "a flat one (no reads) gets haplotypes from its own path like
any other sample" holds for a record without reads, not for a sample without reads anywhere in the
window. Such a sample scores worst in the haplotype choice (no shared records, ½ + 2·√(1/12)), but in
a cohort of 17 samples or fewer every other sample is chosen, so its haplotypes are read.

Related: §7 says an absent template "has shares `1/A` for both haplotypes". The walk of an absent
sample still produces a re-estimated haplotype there, because its likelihood is read as flat
(`records.rs:174–180`). Say whether the store holds `1/A` at absent rows, or whether `1/A` is
applied when the store is read.

**Fix.** Add a rule: a sample not walked in the window carries its `GP` allele shares as both
haplotypes in every round — the prototype's starting value (`phaser.py:263–265`) — and is offered as
a candidate like any other. State the absent-row rule once, at the store.

### I6. The interfaces the coder needs are not stated

**Section:** §3.

- **`walk`.** §3 calls returning the forward vectors "the one change to `model.rs`". The draw also
  needs `persistence` (private, `model.rs:183`), and the re-estimation needs the blurred shares (in
  the private `SiteTerms`). Returning forward vectors from `walk` would also make the fit's 12 grid
  walks (`fit.rs:112–121`) carry them for nothing. Recommend a second function in `model.rs`,
  `walk_and_draw(sites, templates, parameters, generator) -> (WalkResult, two haplotypes per site)`,
  sharing the forward code with `walk`.
- **`refine_sample` has to split in two.** One part is a per-round solve, returning the sample's
  new haplotypes at every walked record, its posterior at its compared records (for the next
  round's choice) and its posterior, prediction and background-copy probability at core records
  (to accumulate). The other is a final step that applies the no-read rules and the fill check to
  the means and builds `SampleSolution`. Name both and say what each returns.
- **The fit runs in round 1 only, and that must be enforced.** `fit_window` adds the window's
  log-likelihoods to running totals (`fit.rs:122`); calling it again in round 2 would count the
  window twice, against haplotype templates. Later rounds must call `SampleFit::parameters`.
- **The background's position.** `settled_genotype`'s frequencies are sliced at
  `templates.len() * allele_count` (`refine.rs:325`); in later rounds that must be the round's
  template count.
- **A one-round configuration for tests.** §7 requires round 1 to reproduce today's pass, but with
  `R` ≥ 2 the output is never round 1's. Make `R` an internal parameter, and define `R` = 1 as
  "write round 1", so the existing `refine.rs` tests run unchanged as the check.
- **The header.** Besides `rounds=R`, write `haplotype_templates=16`: a reader of the VCF otherwise
  cannot tell which model made it.

### I7. The test plan cannot tell a bug from chance, has one test that will be flaky, and misses the cases the prototype never ran

**Section:** §8.

1. **§8.1's "within a tenth" against seed-to-seed spread.** The Rust random generator cannot
   reproduce numpy's, so a Rust run and a prototype run differ by chance as well as by bugs. I ran
   the prototype on the first cohort at three seeds: 689, 724 and 687 wrong; fills 12, 8 and 11
   wrong (`tmp/phasing_review/seeds.out`). The spread is about 5% — half the tolerance. **Add an
   exact test:** make the draw read its uniforms from a stream that can be injected, implement the
   same generator in `phaser.py` (SplitMix64 is a few lines), and require the Rust and Python
   haplotypes and posteriors to agree to 1e-9 after `R` rounds on a small cohort. Then §8.1 can
   be judged against the measured seed spread.
2. **§8.2's draw test is probabilistic as written.** With `ε` = 0.01 every template's blurred share
   is at least 0.005, so no state has zero weight and "the path must follow it everywhere, with any
   seed" fails for some seeds. Run it with `mismatch` = 0 in `WalkParameters`, every other template
   and the background carrying the other allele, and accept either the autozygous state on that
   template or the pair state on it twice.
3. **§8.2's re-estimation example needs copy two's template stated.** With copy one on a
   homozygous-reference template, a likelihood that allows only the heterozygote, and copy two on
   another homozygous-reference template, `P(0,1)` and `P(1,0)` are both 0.995 × 0.005: copy two gets
   the alternative only half the time. The "subtraction" result needs copy two's template to be
   uniform or to lean to the alternative. Give the hand-computed numbers.
4. **Missing items:**
   - multi-allelic and STR records — the prototype read biallelic SNPs only
     (`phasing_ceiling.py:56`, `:69`), so the general disagreement formula
     `Σ_a h(a)·P(genotype without a)` and the re-estimation at `A` > 2 are untested;
   - absent columns, and a sample with no informative record in the window (I5);
   - small cohorts — 2, 3 and 17 samples, where every sample is a template;
   - leave-the-target-out on a few targets, which §7 promises and §8 omits;
   - `GQ` calibration (wrong calls per `GQ` bin, on and off), since the sample's own reads now enter
     twice and §4.2's floor is a `GQ`;
   - §10 Q3's 3 and 8 rounds, as their own item.

### I8. "Window's first core record" does not identify a window, and no generator is named

**Section:** §5, *Determinism*.

A position repeats on every contig, and a record index restarts at each contig, so either alone
gives two windows the same seed. That is deterministic but correlated. The stage already numbers
sites from the run's first site (`stage.rs:130`), which is unique and depends only on the input.
The crate has no random-number dependency (`Cargo.toml` `[dependencies]`); the code uses
hand-written SplitMix64 or xorshift (e.g. `src/cli/estimate_contamination.rs:1005`).

**Fix.** Seed from the run constant, the first core site's run-wide number, the sample and the
round, mixed by SplitMix64. Specify the categorical draw: one uniform in [0, 1) compared with
cumulative weights in a fixed state order — autozygous states first, then pairs row by row, as
`phaser.py:92–96`. This is also what I7.1's exact test needs.

### I9. The prototype measured phasing on biallelic SNPs only; the spec applies it to every record

**Section:** §1, §2.3.

The table in §1 comes from biallelic SNP records (`phasing_ceiling.py:56`). The pass walks every
PASS record, indels and STR records with many alleles included, and phasing changes their
genotypes too. Say this in §1, so the table is not read as covering them, and make I7.4's
multi-allelic test part of §8.

---

## Minor

- **m1. §6 and §8.4 quote the 32-template figure as the adopted result.** Row 7 of §1 gives 115 for
  the adopted 16 haplotype templates and 113 for 32; §6 (line 244) and §8.4 (line 292) say 113. Use
  115.
- **m2. §2.2's "1.9% for round 5 alone".** A rerun with the prototype's own seed gives 35 of 1,645,
  2.1% (`seeds.out`). Other seeds give 1.7%.
- **m3. The oracle's command line.** `phaser.py` defaults to 4 rounds (`phaser.py:243`); the spec's
  numbers need `--rounds 5`. Give the full command and the cell directories
  (`tmp/founder_grid/…`) in §3 or §8. The raw tables behind §1 are not in the repository.
- **m4. Q2's "lower on five" falls within seed noise.** 689 against 699 on the first cohort
  compares one seed with one seed; three seeds of the 16-template run span 687–724. The ruling
  ("within a tenth") is unaffected; the "lower" wording should go.
- **m5. §1.1: "This is what makes the fills more reliable than Beagle's."** No arm isolates soft
  against hard haplotypes. Beagle also phases from hard calls and leaves the target out. Say
  "consistent with", or add a hard-haplotype arm to the prototype.
- **m6. Pass §6's forward-vector figure** ("as `f32`, 44 MB", line 418) disagrees with the code
  (`f64`, `model.rs:226–227`) and with this spec. Fix it in the pass spec when correcting I4.

## Citations checked

All `file:line` citations in the spec point at what they claim: `refine.rs:189`
(`fits.par_iter_mut()`), `model.rs:202` (`walk`), `model.rs:226–274` (the forward vectors and their
loop), `model.rs:236–270` (the forward step), `records.rs:184` (`template_shares`). The named items
exist: `templates::choose_templates`, `templates::disagreement_bound`, `SiteRecord::background_shares`,
`model::GenotypePairs` (crate-visible), `fit::SampleFit::fit_window`, `budget_sample_records` in
`windows.rs`, `07_decision_separate_tool_and_route.md`, and `tomato_fill_test.sh`. The pass spec
sections cited (§2.1's "a thread of its own", §2.3, §3.3's "never a template", §3.4–§3.6, §4.2,
§4.3, §5, §6) say what this spec attributes to them.
