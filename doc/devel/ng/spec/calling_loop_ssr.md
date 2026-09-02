# ng — the STR calling loop: what stands, what is missing, and how the pieces close

*2026-09-02. No new code yet — this settles what remains. Companion documents:
[`run_ssr_observations.md`](run_ssr_observations.md) (upstream: the routing, the walk, and
the merge kind this loop consumes), [`calling_em_loop.md`](calling_em_loop.md) (the loop
itself, written for both paths), [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)
(selection at a tract, settled on paper),
[`read_likelihoods.md`](read_likelihoods.md) §4 (the STR emission),
[`calling_priors.md`](calling_priors.md) §5 (the STR seed),
[`calling_quality.md`](calling_quality.md) §8 (quality at a tract, deferred there),
[`vcf_output.md`](vcf_output.md) (the record).*

---

## 1. What this is

**ng gets its own STR caller — not a port of production's.** Decided 2026-09-02 (owner):
*"we don't want to just port the production STR caller; we want to create an STR caller
using ng's approach, ng's likelihoods and priors, and a new STR variant loop caller. Many of
these pieces are already built."*

That last clause is the finding this document rests on, and it is stronger than the project's
own written record says. [`calling_quality.md`](calling_quality.md) §8 states *"nothing in ng
can score a tract yet"*, because *"the repeat-tract read-likelihood row and the repeat-tract
candidate path are both unwritten"*. **Half of that is stale.** The read-likelihood row is written,
its per-locus parameter assembly is written, the genotype prior's STR seed is written, and
the frequency loop genotypes a repeat tract end to end in its own tests today
([`a_repeat_tract_is_called_from_its_reads`,
`summarise_condition.rs:7623`](../../../../src/ng/calling/inference/summarise_condition.rs)) —
given a candidate set somebody hands it. The candidate path is the half that is genuinely
unwritten, and around it a driver branch, the record's annotation, and a quality decision.

So this document is an inventory with obligations attached: the SNP/indel loop's chain, the
same chain at a tract, what stands at each stage (with the code), and precisely what is
missing. Its scope ends where a `VcfRecord` for a tract exists.

### 1.1 Goals

- A tract cohort observation entering the driver leaves as a written record or a counted,
  named refusal — through ng's loop, scored by ng's stutter emission under ng's prior.
- Every missing piece is named here with its owner, so "build the STR caller" decomposes
  into work items that already have settled designs wherever one exists.
- The measured baseline is stated so done is checkable: on the GIAB benchmark at 30×/50×,
  ng's recall on ground routed as repeat is exactly 0.000 today, against the production
  caller's 0.855–0.990 and freebayes' 0.818–0.874
  ([the report](../../reports/ng_str_path_losses_2026-09-02.md) §5).
- STR allele discovery is built and its effect measured (§3.5): how many alleles it finds,
  and whether it improves the calling — with the on/off default set by that measurement.

### 1.2 Non-goals, and what this document does not do

- **It does not redesign anything already settled.** Selection's design is
  [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md); the loop's mechanics — passes,
  leave-one-out, stopping, the discovery and slippage rounds — are
  [`calling_em_loop.md`](calling_em_loop.md); the emission is
  [`read_likelihoods.md`](read_likelihoods.md) §4; the seed is
  [`calling_priors.md`](calling_priors.md) §5. Where those documents left a question open it
  stays open here, pointed at.
- **Not bundles.** `LocusEvidence::Ssr` can carry a bundle's evidence by type, but no
  generator mints one and no design says what its candidates would be
  ([`run_ssr_observations.md`](run_ssr_observations.md) §8).
- **Not the pre-pass or its CLI.** The STR fits exist in code
  ([`joint/ssr_fit.rs`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs),
  [`slippage_curve.rs`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs),
  [`stratum_fits.rs`](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)) and
  the parameters file carries their output
  ([`bindings.rs:614-629`](../../../../src/ng/calling/parameters_file/bindings.rs)), but no
  command runs a fit yet — so the first STR calls will be `--defaults` calls, scored under
  the shipped stutter constants with `Defaulted` warrants (§3.4). That is honest and weak,
  and giving the fits a command is **deferred future work** (owner, 2026-09-02) — neither
  this document's nor the plan built from it.

### 1.3 Vocabulary

- **The loop** — [`calling_em_loop.md`](calling_em_loop.md)'s iteration: score every
  sample's genotypes from the read likelihoods and a prior built on cohort allele
  frequencies; re-fit the frequencies from the posteriors; repeat to convergence.
- **The ladder** — the tract's candidate structure: sequences pooled into rungs keyed by
  whole-repeat count, ±1 rescue, per-rung support. Selection's output.
- **Stutter / slippage** — a read reporting one or more whole repeat units more or fewer
  than the DNA it came from, from polymerase slippage during copying. The reason a tract
  cannot be scored by the SNP/indel emission.
- **A stratum** — the (period, reference repeat count) cell the slippage numbers are fitted
  and looked up by.
- **Ground** — the stretches of reference sequence a run analyses, in bases; the run
  report's own word.

---

## 2. The chain, side by side

The SNP/indel path's stations, in the order the driver runs them, and the same station at a
tract. **Built** means shipped in `src/` with tests; **settled** means a spec exists and no
code; **missing** means neither.

| station | SNP/indel path | STR path | state |
|---|---|---|---|
| cohort observation stating its kind | the merge | same object, kind carried | **[`run_ssr_observations.md`](run_ssr_observations.md) §4** — upstream of this document |
| candidate selection | [`select_generic`](../../../../src/ng/calling/allele_candidates/generic.rs) | `select_ssr` — the ladder, nomination, ±1 rescue, periodicity verdict | **settled, code missing** (§3.1) |
| evidence shaping | [`shape_generic_locus`](../../../../src/ng/calling/evidence_shaping.rs) | [`shape_ssr_locus`, `evidence_shaping.rs:443`](../../../../src/ng/calling/evidence_shaping.rs) | **built** |
| read likelihood | base-quality emission ([`read_likelihoods.md`](read_likelihoods.md) §3) | stutter + substitution emission, Model A ([`StutterSubstitutionEmission`, `ssr_emission.rs:339`](../../../../src/ng/calling/likelihood/ssr_emission.rs)) | **built** |
| per-locus parameter assembly | read-group calibrations | [`repeat_tract_parameters.rs`](../../../../src/ng/calling/inference/repeat_tract_parameters.rs): a scoring context per (read group, candidate), stated constants where a fit is absent, warrants throughout | **built** |
| genotype prior | fitted spectrum, projected ([`calling_priors.md`](calling_priors.md) §4) | the stratum's length spectrum, indexed from the reference tract length ([`calling_priors.md`](calling_priors.md) §5; `fill_seed_share_per_candidate`) | **built** |
| the frequency loop | [`SummariseConditionLoop`, `summarise_condition.rs:2066`](../../../../src/ng/calling/inference/summarise_condition.rs) | **the same loop** — it branches on the evidence kind internally, and its tract arms are exercised end to end in tests | **built** |
| discovery round | specified, off by default | same design ([`calling_em_loop.md`](calling_em_loop.md) §4.1); constants shipped and inert ([`inference/mod.rs:149-197`](../../../../src/ng/calling/inference/mod.rs)) | **settled, code missing — in scope here, with its measurement (§3.5)** |
| genotype quality | posterior of the called genotype | identical — a property of a posterior, not of a kind | **built** |
| site quality | fold + artifact correction | fold runs; artifact tests deliberately skip a tract ([`summarise_condition.rs:1913`](../../../../src/ng/calling/inference/summarise_condition.rs); test at `:5719`) | **built as far as designed; the design gap is §3.3** |
| the driver | [`call_one_generic_locus`, `callers.rs:813`](../../../../src/ng/run/callers.rs), unconditional | a branch on the observation's kind | **missing** (§3.2) |
| the record | [`assemble_record`](../../../../src/ng/vcf/assemble.rs) | same function; `repeat_tract: Option<TractAnnotation>` already on its input, the run passes `None` ([`records.rs:263`](../../../../src/ng/run/records.rs)); `STR`/`RU`/`PERIOD`/`REPCN` and the tract FILTERs already encode ([`encode.rs:257`](../../../../src/ng/vcf/encode.rs), [`FilterVerdict`, `vcf/mod.rs:684`](../../../../src/ng/vcf/mod.rs)) | **built; two inputs unwired** (§3.2) |

One fact worth a line on its own, found only by opening the file: **every run today already
instantiates the STR machinery.** The genotyper the command builds is
`SummariseConditionLoop<StutterSubstitutionEmission, MarginalizedDirichletPrior>`
([`call_from_alignments.rs:627`](../../../../src/pop_var_caller_exp/call_from_alignments.rs)),
so the STR emission, its scratch, and the loop's tract arms are compiled into every calling
run and reachable the moment an `Ssr`-shaped [`LocusEvidence`,
`calling/mod.rs:392`](../../../../src/ng/calling/mod.rs) arrives. Nothing constructs one.

---

## 3. What is missing, precisely

### 3.1 Selection: `src/ng/calling/allele_candidates/ssr.rs`

The one missing piece of statistics, and its design is done:
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) settles the ladder, which rungs a
sample nominates, the two-spellings rule, the gates deliberately not ported, and the
periodicity verdict; its implementation plan
([`impl_plan/candidate_alleles_ssr.md`](../impl_plan/candidate_alleles_ssr.md)) orders the
build, with the merge's kind field as Milestone A — now owned by
[`run_ssr_observations.md`](run_ssr_observations.md) §4.

What the rest of the chain contracts from it, all three already typed at the consumer:

- a `CandidateAlleles` whose kind is `LocusKind::Ssr(_)` — the parameter assembly asserts
  exactly this
  ([`repeat_tract_parameters.rs:838-847`](../../../../src/ng/calling/inference/repeat_tract_parameters.rs));
- `candidate_repeat_counts`, one `NonZeroU32` per candidate — the field on
  `LocusEvidence::Ssr` whose documentation names its producer as *"repeat-tract candidate
  selection, which is unwritten"*
  ([`calling/mod.rs:416-431`](../../../../src/ng/calling/mod.rs)); not derivable from the
  bases, because an interrupted tract holds fewer whole repeats than its length says;
- the tract FILTER verdicts that are selection's to mint — `NotPeriodic`,
  `TooManyAlleles`, `LowDepth` — for which the record's vocabulary already exists and is
  already declared in every run's VCF header.

### 3.2 The driver branch, and the record's two unwired inputs

`call_one_generic_locus` becomes a dispatch on `CohortObservation::kind` — the generic arm
unchanged, and a tract arm that runs: `select_ssr` → `shape_ssr_locus` → the same
`genotyper.call_locus` → `assemble_record`. The scratch is already in place (the run's
`CallingScratch` is parameterised by the STR emission's own scratch type,
[`callers.rs:551`](../../../../src/ng/run/callers.rs)). Two values the record's assembly
already accepts must start arriving:

- `repeat_tract: Some(TractAnnotation::new(motif))` from the observation's kind — the field
  and its `None` placeholder are at
  [`assemble.rs:79-84`](../../../../src/ng/vcf/assemble.rs) and
  [`records.rs:263`](../../../../src/ng/run/records.rs);
- the FILTER verdict from selection/loop rather than the converged-or-not derivation alone.

Same station, one more obligation: **the run report.** A tract called, a tract refused
`TooManyAlleles`, and a tract on unroutable ground are three different facts, and the report
that today prints *repeat tracts this caller has not built yet* must partition them the way
it already partitions the generic path's outcomes.

### 3.3 Quality at a tract: the risk is stated, and the decision waits for a measurement

What exists: the genotype quality is kind-blind and runs; the site-quality fold runs at a
tract; the strand/read-position artifact tests skip it on purpose — slippage, not strand, is
what goes wrong at a tract, and slippage is already inside the emission
([`calling_quality.md`](calling_quality.md) §8's reasoning, which survives even though its
"nothing can score a tract" premise did not).

**The possible problem, stated so nobody rediscovers it in a benchmark.** A tract's QUAL is
the fold's claim that the cohort is not all homozygous-reference, and at a tract that claim
leans entirely on the stutter model pricing slip products correctly. If the model under-prices
them anywhere — a stratum scored from the shipped defaults rather than a fit (§3.4), a
homopolymer at high depth where many independent reads all slip the same way — the fold does
at tracts what the uncorrected SNP baseline did at depth: **grows confidently wrong**, because
every extra read compounds the same mispriced error. The SNP path went through exactly this —
its QUAL inflated with depth until two artifact tests were bolted on, and the bolting-on
happened only after the measurement
([`qual_fp_depth_inflation_2026-06-10.md`](../../reports/qual_fp_depth_inflation_2026-06-10.md))
— and production's STR caller ended somewhere stronger still: it does not trust a posterior
for the emission decision at all, and its bake-off picked a port of freebayes' formula over
its heuristic and its likelihood-ratio arms
([`ssr_emission_model_comparison_2026-07-08.md`](../../reports/ssr_emission_model_comparison_2026-07-08.md)).
Neither history proves ng's fold fails at tracts — ng's stutter emission is a different model
than either measured — but both say the naive answer has failed here before, in both
directions worth watching: confident false positives at depth, and a QUAL too timid to gate
on at three reads.

**The options, as arms of one experiment:**

- **A — the inherited fold as it stands.** What ships first.
- **B — the fold plus a tract-specific correction** in the slot the artifact correction
  occupies on the SNP path, its test designed against whatever failure shape arm A actually
  shows. Deliberately not designed here: designing the correction before the miscalibration
  is measured is how an invented test gets built.
- **C — a per-locus emission decision instead of the posterior**, production's shape.
  Production's own caller is the comparator for this arm — run on the same ground, not
  ported; `src/ssr/` stays frozen.

**The obligation — decided 2026-09-02 (owner): the tract QUAL design is chosen only after an
experiment observes the options' behaviour.** The experiment, so the plan can carry it as a
milestone rather than a sentiment:

- **Data:** the GIAB tract ground at 30× and 50× (the split and scripts of
  [the loss report](../../reports/ng_str_path_losses_2026-09-02.md) §5), with production and
  HipSTR beside ng from the standing benchmarks; and the STR simulator for exact truth at
  settable slippage, which the bake-offs harness already uses.
- **Measured, per arm:** calibration — records binned by QUAL against the share truly
  variant, so "QUAL 30" can be checked against "wrong about one time in a thousand"; and
  gateability — precision and recall as a QUAL threshold sweeps, the instrument the
  qual-analysis work already built for SNPs. Split homopolymers from period-2+ tracts, and
  fitted-parameter cells from `Defaulted` ones, because those are the axes the risk names.
- **The decision rule:** arm A ships as the tract QUAL if it is calibrated and gateable to
  the same standard the corrected SNP QUAL reaches on the same benchmark; where it is not,
  the measured failure shape is what `calling_quality_ssr.md` designs arm B against, with
  arm C's numbers as the bar it has to beat.

**Home for the final design:** `calling_quality_ssr.md`, exactly where
[`calling_quality.md`](calling_quality.md) §8 already put it — now written *from* the
experiment's report rather than before it. Until then the inherited site quality ships,
labelled by this section as unvalidated at tracts. **The STR loop's implementation plan must
carry this experiment as a step** — it is registered meanwhile in
[`calling_bakeoffs.md`](../impl_plan/calling_bakeoffs.md)'s next-plans list so it cannot be
dropped.

### 3.4 Parameters: what a tract is scored from, until the fit work happens

**Not a missing piece of this work — context for reading its first results.** A tract is
scored from per-stratum slippage, a length spectrum, and substitution rates; where a
parameters file supplies fitted values the assembly uses them, and every absence is answered
with a stated constant and a warrant — the shipped HipSTR stutter model, `Defaulted`, per
cell
([`repeat_tract_parameters.rs:22-46`](../../../../src/ng/calling/inference/repeat_tract_parameters.rs))
— with the outlier weight production's 0.01, declared inherited. Producing fitted files is
**deferred future work** (§1.2, §6), so the first tract calls will be `--defaults` calls:
uncalibrated slippage the way the path already calls SNPs with uncalibrated base qualities —
honestly labelled by the warrants and the run report, and measurably weaker. Every
acceptance number in §8 is read with that in mind.

### 3.5 Allele discovery is in scope — build it, and measure what it buys

**Decided 2026-09-02 (owner): implementing STR allele discovery, and measuring its effect,
is part of this work.** The two questions the measurement answers are the owner's own: *do
we discover many alleles?* and *does it improve the STR calling?*

**What discovery is.** The candidate set a tract is called over comes from what the reads
showed — but a true allele can hide *under* stutter: every read carrying it is booked as a
slip product of a called length, so its repeat count never surfaces as a candidate. The
discovery round looks after convergence for lengths the posteriors keep explaining away,
admits the ones that clear an evidence bar, and reruns the loop to a fixpoint — HipSTR's
mechanism, set out arm by arm in [`calling_em_loop.md`](calling_em_loop.md) §4.1. That
design is settled, its evidence-bar constants are shipped and inert
([`inference/mod.rs:158`](../../../../src/ng/calling/inference/mod.rs) — *"inherited from
HipSTR's high-depth human setting and soft"*), and on the STR path an admitted allele changes
every genotype's likelihood — the rows are recomputed rather than carried over
([`calling_em_loop.md`](calling_em_loop.md) §4's table), while the per-read emission columns
append. The round's body is what is unwritten.

**Sequencing.** Discovery grows a candidate set, so it runs only where a candidate set
exists: after §3.1's selection and §3.2's driver branch. It is this scope's last step, not
its first.

**The measurement, so "does it help" is a report and not an impression** — the shape is
Q3's, already specified as the bake-offs plan's F3 report, run here on the STR path:

- **How many alleles.** How often the round fires at all (loci in ten thousand — if a
  handful, everything after is a rounding error); on firing loci, alleles admitted, alleles
  surviving the prune, rounds to fixpoint.
- **Does it improve the calling.** Discovery off against on, same ground, same truth:
  recall and genotype concordance on GIAB tract ground at 30×/50×, exact-truth accuracy on
  the STR simulator at settable slippage — and the danger case, tomato at three reads,
  where the evidence bar's read floor binds differently and a stray read minting an allele
  is the failure to watch (the bar sweep F3 specifies).
- **The default is the measurement's to set.** Discovery ships built but off; it becomes
  the default only if the on-arm wins the comparison above. Its cost — the genotype
  likelihoods recomputed each round — is reported beside its benefit.

**The per-locus slippage re-fit stays out.** It is the other optional round
([`calling_em_loop.md`](calling_em_loop.md) §5.1), it moves fitted error numbers rather
than the candidate set, and nothing in the owner's ruling touched it — deferred, §6.

---

## 4. One sample and a thousand, three reads and three hundred

- **One sample.** Selection's bar and the periodicity verdict are cohort-denominated but
  defined at N = 1 ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §9); the seed is
  the fitted spectrum, which at `--defaults` is a stated constant — so a single-sample tract
  call leans hardest on the emission, and says so through its warrants. No-calls are the
  honest common case at low depth, and the tract `LowDepth` FILTER exists for the locus-level
  version of the same fact.
- **A thousand samples.** The loop's cost model is unchanged — the tract arms feed the same
  frequency loop, and the parameter table is (read groups × candidates), not
  (samples × anything). The allele cap binds harder at a tract (many rungs segregate at
  mutation-hot loci); `TooManyAlleles` is the named refusal, and cutting the worst-evidenced
  rungs while keeping the locus is the same policy the generic path already applies.
- **Three reads / three hundred.** The stutter emission's shape is per read; depth changes
  confidence, not the model. The read cap (1,000/locus, soft) is the only depth knob and
  sits far above 300.

## 5. Reuse map

| what | source | how |
|---|---|---|
| the loop, both optional rounds' design | [`calling_em_loop.md`](calling_em_loop.md) + [`summarise_condition.rs`](../../../../src/ng/calling/inference/summarise_condition.rs) | as is — already branches on the evidence kind |
| emission Model A | [`likelihood/ssr_emission.rs`](../../../../src/ng/calling/likelihood/ssr_emission.rs), [`ssr.rs`](../../../../src/ng/calling/likelihood/ssr.rs) | as is; why this model is [`read_likelihoods.md`](read_likelihoods.md) §4.1 |
| per-locus parameter assembly | [`repeat_tract_parameters.rs`](../../../../src/ng/calling/inference/repeat_tract_parameters.rs) | as is |
| STR seed | [`calling_priors.md`](calling_priors.md) §5 + `genotype_prior` | as is |
| selection design | [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) (+ arch, + plan) | build `allele_candidates/ssr.rs` to it |
| production's selector | `src/ssr/cohort/candidate_set.rs` (frozen) | **differential oracle only** — with its three replaced rules switched in, `select_ssr` must reproduce it (the plan's own acceptance); never a dependency |
| record + fields | [`vcf/assemble.rs`](../../../../src/ng/vcf/assemble.rs), [`encode.rs`](../../../../src/ng/vcf/encode.rs) | as is; wire the two inputs |
| accuracy harnesses | `benchmarks/giab/` (this year's three-caller split), `benchmarks/ssr_hg002/` (HipSTR/truth eval), `benchmarks/ssr_tomato1/` | acceptance, §8 |

## 6. Deferred, with a recommended home

- **`calling_quality_ssr.md`** — the tract QUAL decision (§3.3), after first records exist.
- **Bundles** — [`run_ssr_observations.md`](run_ssr_observations.md) §8's deferral; a
  bundle's *evidence* type already fits `LocusEvidence::Ssr`, its candidates are the open
  half.
- **The per-locus slippage re-fit** — the optional round that moves the fitted error
  numbers at a locus; stays with
  [`impl_plan/calling_bakeoffs.md`](../impl_plan/calling_bakeoffs.md) (Milestone D, report
  F2). Discovery is no longer deferred — it is §3.5's scope.
- **The fit-mode command** — deferred future work (owner, 2026-09-02): no step of this
  spec or of the plan built from it. Home: the run driver's fit mode, when that work is
  scheduled; §3.4 says how tract calls are read until then.

## 7. Resolved decisions and open questions

- **ng's own caller, not a port — decided** (owner, 2026-09-02). Production's STR caller is
  a comparator and a differential oracle for selection; nothing in `src/ssr/` becomes a
  dependency, per the standing freeze.
- **One loop for both paths — decided long since and standing**: the tract path enters
  `SummariseConditionLoop`, it does not get a second loop
  ([`calling_em_loop.md`](calling_em_loop.md) §2; the code already holds this).
- **OPEN — the tract's QUAL, and it may only be closed by §3.3's experiment** (owner,
  2026-09-02): inherited fold vs a tract-specific correction vs an emission-style decision.
  Leaning: ship the inherited fold, run the experiment on the first real tract records,
  decide in `calling_quality_ssr.md`. Until it closes, tract QUAL is not comparable across
  callers and the file's own labelling says so.
- **OPEN — `max_candidate_alleles` at a tract**: the run-wide cap (default 6, reference
  included) was set on SNP/indel grounds; tracts routinely segregate more rungs. Leaning:
  keep one cap until the `TooManyAlleles` counts from real runs say otherwise; the run
  report will carry the number.

## 8. How we know it works

The baseline is exact zeros, so the first measures are blunt and the later ones are the
project's standing harnesses:

1. **It calls at all**: on GIAB at 30×/50×, recall on repeat-routed ground moves off 0.000.
   The bar to state beside it: production 0.855/0.909 (indels, 30×/50×) and 0.990 (SNPs);
   freebayes 0.818–0.874. Same split, same truth, same script
   (`benchmarks/giab/src/`, [the report](../../reports/ng_str_path_losses_2026-09-02.md) §5).
2. **Selection parity**: the differential against production's selector, as
   [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §13 specifies — failing states at
   both ends.
3. **Genotypes, not just sites**: the GIAB genotype-concordance panel already in the
   dashboard, where the SNP/indel path's known weakness (indel genotypes wrong one time in
   five) sets the number to beat at tracts.
4. **The record**: a tract record round-trips `STR`/`RU`/`PERIOD`/`REPCN` and its FILTER
   through `bcftools` untouched — the encode tests carry most of this today.
4b. **The QUAL experiment ran before the QUAL decision** (§3.3): its report exists, per arm,
   with the calibration and sweep numbers split by period and by parameter provenance —
   done is a written report, not a chosen answer.
4c. **Discovery's report exists** (§3.5): fire rate, alleles admitted and surviving, and the
   off-against-on accuracy comparison on GIAB, the simulator, and tomato at three reads —
   and the shipped default matches what the report says.
5. **Cohort invariance**: E2's byte-identity across thread counts, with tract records in
   the file.
