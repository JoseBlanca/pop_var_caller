# ng — synthetic validation (the independent tester)

*Status: design spec (2026-07-29). **No code yet — this settles the design**, and it is
written to be argued with. Companion arch and plan docs do not exist yet. Grounded in a
source read of `src/ssr/cohort/sim.rs`, `src/ng/locus_generation/mod.rs`,
`examples/ssr_delimiter_comparison.rs`, `examples/ng_synthesize_stress_reads.rs`, and the
ng blocks of `PROJECT_STATUS.md`.*

*Naming: **STR** in prose, `ssr` in code.*

---

## Why this exists

ng's tests today answer **"does this code do what its author meant?"** They are good at it —
the parity harnesses, the leftmost-property oracle, the 200,000-case soaks. What none of them
answers is **"does this code do what the feature requires?"**, because almost every one of
them is written by the same hand, from the same reading of the spec, as the code it checks.

A variant caller has a way out that most software does not: **truth is constructible.** If we
synthesize reads from a genotype we chose, we know the answer before we run the caller. An
independent tester that builds inputs this way can check the *feature*, not the
implementation's self-consistency.

The design question this doc settles is not "should we simulate data" — the repo already does,
in three places. It is **what makes a simulator an oracle instead of a mirror**, and that turns
on one rule (§2) that the existing simulator deliberately does not follow.

**Jargon, once.** *Constructed truth* is the genotype/allele we injected, known before the run.
A *metamorphic relation* is a property linking two runs — transform the input a known way, and
the output must change a known way — which needs no truth at all. *Calibration* is whether a
reported confidence matches the observed error rate. *Defect injection* is deliberately
breaking the code to prove a test can fail. A *shrunk case* is a failing input a property-test
framework has minimised to the smallest form that still fails.

---

## 1. Goals, non-goals, and what this is not

**Goals.**

1. **Constructed truth, independently generated.** Inputs whose correct output is known by
   construction, produced by code that does not share a model with the code under test (§2).
2. **Systematic feature coverage.** A declared axis matrix (§4), swept and *logged*, so
   "we covered indels" is a fact rather than a vibe.
3. **Failures arrive as reproducers.** A failing case is a seed plus a shrunk input plus a
   runnable command — never prose. Triage is running a test, not reading an argument.
4. **Every harness is proven able to fail** (§5), by injecting the defect it claims to catch.

**Non-goals** — reasonable to want, deliberately excluded:

- **Replacing ng's existing tests.** The parity harnesses check byte-identity to production;
  this checks agreement with the *feature*. They fail on different things and both are wanted.
- **Performance.** No timing assertions, no benches. `benches/ng_ssr_delimiter_perf.rs` owns that.
- **Testing production (`src/` outside `src/ng/`).** Production is frozen — ng is a
  from-scratch caller and steps port back only after the experiments (owner, 2026-07-16).
  Several metamorphic relations in §3.2 would apply to it unchanged; deferred, with a home (§10).
- **Being the calibration study.** Calibration numbers are a research report's output. This
  harness may *compute* them (§3.3), but a calibration drift is a finding for a human, not a
  red test.

**What it does not do.**

- It does not judge whether a *model* is right for tomato — only whether the code implements
  the model the spec describes. A wrong-but-faithfully-implemented model passes, by design;
  catching that needs real data (§6).
- It does not read the implementation to decide what "expected" means. Expected comes from the
  spec, or the test is circular.

---

## 2. The independence rule — the one decision everything else rests on

**A generator that imports the model it is testing cannot detect a wrong model.** It can only
detect a wrong *implementation* of that model, and it will report perfect agreement in exactly
the case that matters most.

This is not hypothetical. The repo's existing cohort simulator says so in its own header:

> *"§1 — the generative side of the same model the kernel scores"*
> — [src/ssr/cohort/sim.rs:4](../../../../src/ssr/cohort/sim.rs)

and it imports `MAX_SLIP`, `PerBaseError`, `StutterLevel`, `StutterShape` from
`param_estimation` — the production parameter types
([sim.rs:36](../../../../src/ssr/cohort/sim.rs)). **That is correct for what sim.rs is for**:
it exists so the EM's "recover what we put in" milestones have a known input, and for that
job sharing the parameterization is the point. It is not an independent oracle and does not
claim to be.

The same shape has already produced a recorded near-miss in ng: the 2026-07-24 normalizer
screen compared ng's three normalizers **to each other**, never against production's raw-byte
behaviour — so it could not have found the soft-mask defect that a later, differently-framed
comparison did (PROJECT_STATUS, normalization plan + step-2 blocks).

**The rule.** The generator for a given feature **may not import from the module under test**.
Concretely, for the STR observation path it may not use:

| forbidden to the generator | why |
|---|---|
| `ng::alignment::StutterModel`, `StutterRates`, `MAX_SLIP` | the slip model is what we are checking |
| `ng::alignment::PerQualityEmission`, `FlatEmission` | the error model is what we are checking |
| any `BestPathAligner` / `MarginalAligner` impl | the thing under test |
| `ng::region_typing` tract delimitation | truth is where we *put* the tract, not where the scanner finds it |

It **may** use: coordinate and identity types (`ContigId`, `GenomeRegion`, `Motif`, `Bp`), the
BAM/FASTA writers, and `RefSeq` to read back a reference it wrote. Those carry no model.

**The price, stated plainly.** The generator has to express stutter and sequencing error a
second time, independently. That is duplicated work and it will drift from the production
model — **and the drift is the signal**. Two independent expressions of the same biology
agreeing is evidence; one expression checked against itself is not.

**Soft, and worth arguing about (§11, Q2):** the strongest form of this rule is that the
model-bearing parts of the generator are written by someone — or some agent — who has read the
spec and *not* the implementation. That is cheap to arrange and expensive to verify.

---

## 3. What gets asserted — three families, kept apart

They are separated because they have different failure semantics, and mixing them produces a
suite that flakes and then gets muted.

### 3.1 Constructed truth — exact, per-case, hard

Inject a known truth, run the code, assert recovery. **These are exact assertions**: at
sufficient depth and quality, the recovered allele either equals the injected one or the code
is wrong. No tolerances.

The assertion surface today is `SampleLocusObservations`
([src/ng/locus_generation/mod.rs:34](../../../../src/ng/locus_generation/mod.rs)) — see §8 on
why that is the end of the line, and what changes when a genotyper lands. Per locus:

- the injected allele appears in `observed_sequences`, with support equal to the number of
  reads generated from it (modulo reads the cap discarded — `reads_discarded_by_cap` is
  non-zero exactly when that happened, so the harness can tell the two apart);
- no allele appears that was never injected and is not explained by an injected error;
- `reference_bases` matches the reference the generator wrote;
- `reads_without_observation` accounts for every read that was generated but not observed.

**The negative case is the important one and it is the one ng cannot currently reach.** Inject
*nothing* — pure reference reads with realistic error — and assert nothing variant-like
survives. That is the false-positive axis, and it is production's known-weak dimension:
false-positive SNP QUAL climbs with coverage where freebayes' stays flat
([qual_fp_depth_inflation_2026-06-10.md](../../reports/qual_fp_depth_inflation_2026-06-10.md)).
At the observation seam
"nothing is called" is not expressible, because observations are not calls. §8.

### 3.2 Metamorphic relations — no truth needed, and they work on real data

These link two runs of the *same* code and need no generative model at all, so they sidestep §2
entirely — and they run on real CRAMs as happily as on synthetic ones. **Per unit of effort
this is the cheapest family, and it is the one to build first** (§11, Q1 is where it lives,
not whether).

The set below is what ng's current surface admits. Each is a hard, exact assertion.

| relation | transform | expected |
|---|---|---|
| **translation** | shift a locus and its reads by *N* bases | identical observations, coordinates + *N* |
| **reverse complement** | RC the reference and all reads | identical observations, RC'd and mirrored |
| **sample permutation** | reorder samples in the cohort walk | per-sample observations byte-identical |
| **sample duplication** | run the same sample twice under two names | two identical observation sets |
| **read shuffling** | permute read order within a locus | identical observations |
| **region partition** | split a region set in two, run each, concatenate | identical to one run over the union |
| **reference-only addition** | add reads matching the reference exactly | no *new* non-reference observation appears |
| **depth subsample** | drop half the reads | support counts fall; no new allele appears |
| **soft-mask** | lowercase the reference | identical observations |
| **determinism** | run twice, same input | byte-identical output |

Two of these have already bitten this codebase, which is the argument for the list rather than
for any one entry: the staged-producer multi-interval bug dropped nearly all calls on a
multi-interval single-sample run (*region partition*), and production's raw-byte left-alignment
does nothing at all on a soft-masked reference (*soft-mask*) — a defect ng fixes and that the
step-2 parity fixture demonstrated (PROJECT_STATUS, step-2 block).

**A relation that fails names its own bug.** "Translation-invariance broke at offset 4,096" is
a coordinate-arithmetic bug with its own reproducer attached — there is nothing to triage.

### 3.3 Calibration — aggregate, banded, and never a red test on its own

Over ≥10⁴ sites: does a reported confidence match the observed error rate? Does a point
estimate carry bias? `examples/ssr_delimiter_comparison.rs` already establishes the right
vocabulary for a *ruler* rather than a probability — accuracy, signed bias, spread
([ssr_delimiter_comparison.rs:9-25](../../../../examples/ssr_delimiter_comparison.rs)).

**These never assert per case.** They emit a table; a drift is a finding for a human. Wiring a
tolerance band into CI here is how a suite earns a reputation for flaking, and a suite with
that reputation is off.

---

## 4. Variability — the axis matrix, and where mutation belongs

"Add mutation and variability" has a precise form: **property-based testing over a domain
generator.** `proptest` is already a dependency ([Cargo.toml:126](../../../../Cargo.toml)) and
in use in 13 modules.

**The reason to use proptest rather than a hand-rolled loop is shrinking**, and it is the whole
reason. A failure at 50 samples × 60× depth × a 40 bp tract is not actionable; the same failure
shrunk to 3 reads and 1 sample is a morning's work. That difference decides whether a finding
gets fixed or gets filed.

### The axes

Random generation systematically under-samples the corners, and this repo's record is that the
corners are where the bugs are: B1's zero-flank case ("the gap was in the fixture, not the
code"), C1's band failure at case 1745. So the axes are **declared and swept**, with randomness
*within* a cell rather than instead of the matrix.

| axis | values | note |
|---|---|---|
| marker | STR · generic SNP · indel | generic path is thin until the pileup lands (§8) |
| period | 1, 2, 3, 4, 5, 6 | measured on the tomato cohort: stutter is a mono/di phenomenon, tri–hexa flat to 30 bp |
| tract length | 5 … 60 bp | stutter onset ~10 bp mono, ~15 bp di |
| tract purity | pure · 1 interruption · 2+ | interruptions are a recorded recall gap |
| flank available | both · left-only · right-only · neither | the four `RepeatSpan` cases |
| depth | 1, 3, 10, 30, 100, cap+1 | 3 is the tomato cohort's real depth |
| genotype | hom-ref · het · hom-alt · multi-allelic | het is where the lone-carrier tension lives |
| error rate ε | 0, 0.001, 0.01, 0.05 | 0 is the "must be exact" tier |
| read length vs tract | spanning · partial · flanking-only | drives `ReadCoverage` |
| reference edge | interior · contig start · contig end | off-by-one territory |
| soft-masking | none · tract · flank | tomato SL4.0 has ~227,170 lowercase bases |
| samples | 1 · 2 · 50 | 50 only in the slow tier |

**Soft, all of it.** These values are starting points chosen to straddle known thresholds, not
measured optima. They are the most movable thing in this document.

**Sweep logging is mandatory.** A run records which cells it swept and which it skipped. A
harness that silently caps its own coverage reads as "we covered everything" when it did not
— the same failure `PROJECT_STATUS` flags as "no silent caps" in the review skill's vocabulary.

### Determinism

Every case is a pure function of `(seed, axis-cell, index)`. Both existing generators use a
dependency-free SplitMix64 for exactly this reason — "so the whole comparison is reproducible
from its seed" ([ssr_delimiter_comparison.rs:34](../../../../examples/ssr_delimiter_comparison.rs))
and "so the simulator's determinism does not rest on an external crate's stability"
([sim.rs:46](../../../../src/ssr/cohort/sim.rs)). Reuse that shape; a failure that cannot be
replayed from its seed is not a finding.

---

## 5. The harness must be proven able to fail

**This is not a nice-to-have; it is the single most load-bearing rule in the document**, and it
is here because of this repo's own record. Every one of these is from `PROJECT_STATUS.md`:

- alignment B1 — "every one of B1's 3 Blockers and B3's Majors was a test that could not fail
  — found by mutating the source, not by reading it";
- A3 stutter model — the fixtures "gave paired parameters equal values, so two transpositions
  (`in_up`↔`in_down`, `in_geom`↔`out_geom`) and **the entire out-of-frame direction split**
  were invisible to all twelve tests";
- normalization A2 — "1 Blocker + 4 Majors, every one a test that could not fail";
- C1 banding — every band term "proven load-bearing" by mutation, the failure landing at a
  named case (1745, 288).

A synthetic harness is *more* exposed to this than a unit test, not less: it is large, it is
green by default, and its green is very reassuring.

**The rule.** Each axis ships with at least one named defect it is claimed to catch, injected
and confirmed to turn the harness red. Recorded as a table in the harness's own report:

| axis | injected defect | expected to fail |
|---|---|---|
| period | transpose `in_up` / `in_down` | stutter direction asymmetry at period 1 |
| flank available | swap left/right flank lengths | `RepeatSpan::FromLeft` vs `FromRight` cases |
| reference edge | off-by-one on the window start | contig-start cells only |
| tract purity | route an interruption as in-frame stutter | interrupted cells only |
| translation (§3.2) | drop the offset from a coordinate | all cells |

The precedent is already in the tree: `ng_synthesize_stress_reads.rs` exists specifically "to
show `ng_normalizer_screen` is genuinely discriminating"
([ng_synthesize_stress_reads.rs:1-3](../../../../examples/ng_synthesize_stress_reads.rs)), and
it produced 28/64 disagreements against a real-data screen's 0. Generalize that from one screen
to every axis.

---

## 6. Real data is the validity anchor

Synthetic data sweeps the space. It cannot tell you the space is the right shape — real reads
carry error modes nobody models: optical duplicates, adapter read-through, reference bias,
mapping artifacts at paralogs.

The anchor already exists: `benchmarks/ssr_hg002/` holds the GIAB HG002 TR benchmark v1.0.1 —
our `.cat` catalog, the 300× Illumina BAM, the assembly-truth VCF, and the Tier BED. Tomato
gives the low-depth cohort shape that HG002 cannot.

**Division of labour:**

- **synthetic** — coverage and corners; the only way to reach `depth = 1`, `ε = 0`, or a
  contig-start tract on demand;
- **real** — validity; the only way to find out that an axis is missing.

**Both directions of disagreement are informative.** Synthetic-only failure is a generator
artifact until proven otherwise. Real-only failure means the generator is missing an axis —
and that is itself a finding worth filing, because it names a phenomenon the model does not
represent.

---

## 7. Cross-cutting concerns

**Runtime.** The measured shape of ng's walk: the delimiter is 51.1% of a tomato run, noodles
CRAM slice-MD5 13.0%, block decode 5.9%, region typing 5.2%; the cohort stutter probe hit a
~38 ms/locus/sample wall (PROJECT_STATUS). So the harness needs **two tiers**: a fast tier that
constructs reads in memory and calls the delimiter directly (no BAM, no CRAM, no region typing
— that is the shape `ssr_delimiter_comparison.rs` already uses), and a slow tier that writes a
real reference and a real indexed BAM and runs the whole spine (the shape
`ng_synthesize_stress_reads.rs` uses). The fast tier is the inner loop; the slow tier runs on
demand and covers exactly the code the fast tier bypasses — I/O, region typing, coordinate
plumbing, which is where several §3.2 relations live.

**Errors.** A generator failure must be loud. A case that silently fails to generate and is
skipped reads as a pass, which is the §5 failure mode wearing a different coat. Generation
errors abort the run.

**Concurrency.** None. ng's walk is single-threaded per worker and the harness inherits that.
Cases are independent, so parallelism is available later at zero design cost.

**Memory.** The fast tier holds one locus at a time. The slow tier's footprint is the walk's,
already characterized.

---

## 8. Scope today: ng has no end, and that shapes the harness

**The load-bearing fact for scheduling.** ng's built spine stops at
`SampleLocusObservations` — there is no genotyper, no likelihood over genotypes, and no VCF
writer in `src/ng/`. The module list is `alignment`, `locus_generation`, `read`,
`region_typing`, `ref_seq`, `reference_info`, `tandem_repeat`, `types`.

So the natural end-to-end assertion — *inject genotype G, assert the VCF says G* — **is not
available**, and will not be until the genotyping steps land. What is available today:

| assertion | available now? |
|---|---|
| injected allele appears in `observed_sequences` with the right support (§3.1) | ✅ |
| every metamorphic relation in §3.2 | ✅ |
| delimiter measures the injected tract length (fast tier) | ✅ — this is `ssr_delimiter_comparison.rs`, generalized |
| read filtering keeps/drops the reads we designed it to | ✅ |
| **injected genotype recovered** | ❌ no genotyper |
| **injected nothing → nothing called** (the FP axis) | ❌ observations are not calls |
| QUAL calibration | ❌ no QUAL |

**Design consequence, and it is the one thing here that constrains the coder:** the truth table
is written to carry the genotype from the start, even though nothing reads it yet. Truth is
*"sample S has alleles (A₁, A₂) at locus L, from which these reads were drawn"* — the
observation-level assertions project out of that. Store only observations now and the whole
generator gets rewritten when the genotyper lands.

The two ❌ rows are the most valuable assertions in the document. They are also the reason this
harness's ambition should be re-read against the ng roadmap before a large investment: **half of
what it is for cannot be built yet.**

---

## 9. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| deterministic PRNG | `SplitMix64` in [sim.rs:46](../../../../src/ssr/cohort/sim.rs) and [ssr_delimiter_comparison.rs:34](../../../../examples/ssr_delimiter_comparison.rs) | copy the shape (dependency-free determinism), not the file |
| synthetic reference + indexed BAM writing | [ng_synthesize_stress_reads.rs](../../../../examples/ng_synthesize_stress_reads.rs) | the slow tier's IO scaffolding, lifted out of the example |
| fast-tier read construction + truth scoring | [ssr_delimiter_comparison.rs](../../../../examples/ssr_delimiter_comparison.rs) | the accuracy/bias/spread vocabulary and the loop shape |
| in-memory reference | `InMemoryRefSeq` ([src/ng/ref_seq.rs](../../../../src/ng/ref_seq.rs)) | the fast tier's reference, no FASTA on disk |
| shrinking + case generation | `proptest` ([Cargo.toml:126](../../../../Cargo.toml)) | the `Strategy` for a cell of the axis matrix |
| the assertion target | `SampleLocusObservations` ([locus_generation/mod.rs:34](../../../../src/ng/locus_generation/mod.rs)) | asserted against, not reused |

**No parity oracle is named, because this is not a port.** Its correctness rests on §5
(defect injection) instead — which is why §5 is not optional.

The existing `delimit_parity.rs` and `leftmost_property.rs` are **siblings, not inputs**:
parity checks byte-identity to production, the property oracle checks a definition, and this
checks agreement with constructed truth. Three different questions.

---

## 10. Deferred, with a recommended home

- **Metamorphic relations applied to production** (`src/pileup`, `src/var_calling`) — several
  in §3.2 apply unchanged and production is where the users are. Home: a production-side
  harness after ng's is proven, or the port-back of whichever ng step lands first. Not here:
  production is frozen (§1).
- **Cross-caller differential** (ours vs HipSTR/GangSTR/freebayes on synthetic truth) — a
  research report's job, not a test's. Home: `doc/devel/reports/research/`; the tooling exists
  in `benchmarks/ssr_hg002/scripts/`.
- **Calibration bands in CI** — deferred until there is a QUAL to calibrate (§8) *and* a
  measured baseline to band around. Home: a follow-up spec once the genotyper lands.
- **Continuous/agent-driven operation** — running this harness in a loop that reports
  divergences without a human present. It is the natural consumer, but it is a separate
  concern with its own failure modes (deduplication, budget, ledger). Home: its own spec.

---

## 11. Open questions

Every one of these is genuinely open, with a leaning; **confirm before code.**

**Q1 — Where does it live?** Options: `src/ng/validation/` as a `#[cfg(test)] pub(crate)`
module (reachable from unit tests, invisible to release); `tests/` (integration, but then the
fast tier cannot reach crate internals); `examples/` (where both existing generators live, but
examples are run by hand and rot). **Leaning: `src/ng/validation/` for the generator and the
fast tier, plus one thin `examples/` driver for the slow tier** — matching where
`delimit_parity.rs` and `leftmost_property.rs` already sit. Note the cost, which C1 already
paid: `src/ng/read/input/test_fixtures.rs` is `#[cfg(test)] pub(crate)`, and that is precisely
why no bench can reach the synthetic indexed-BAM builder (PROJECT_STATUS, M1).

**Q2 — How hard is the independence rule?** §2 forbids importing the model. The stronger form
— the generator's model-bearing code is written by an author who has read the spec and never
the implementation — is cheap to arrange and hard to verify. **Leaning: enforce the import ban
mechanically (it is a grep), and treat authorial independence as a strong preference for the
stutter and error models specifically.**

**Q3 — Which family is built first?** **Leaning: metamorphic (§3.2)**, because it needs no
generator, no model, and no truth table, runs on data that already exists, and two of its ten
relations map onto bugs this codebase has already had. Constructed truth (§3.1) is the bigger
prize and the bigger build.

**Q4 — What settles "the generator is good enough"?** The honest answer is §5's injection table
plus §6's real-data cross-check, and neither is a number. **Leaning: the exit criterion is the
defect table, not a coverage percentage** — a harness that catches every injected defect on
every declared axis is done, and one that catches 90% names the gap.

**Q5 — Do we generate reads, or alignments?** Generating *reads* and aligning them with bwa
would test the whole chain including the mapper, but the mapper is not ours and its output is
not reproducible across versions. Generating *alignments* directly (CIGAR + sequence, as both
existing generators do) tests only our code. **Leaning: alignments**, with the mapper's
contribution left to real data (§6) — but this is the axis on which the harness is most
obviously blind, and it deserves an argument.

**Q6 — How much does §8 change the plan?** Half the value is behind a genotyper that does not
exist. **Leaning: build the fast tier and the metamorphic relations now (both are useful
today), design the truth table for genotypes (§8), and defer the constructed-truth cohort
generator until there is something to genotype.** The manager's version of this question is
whether the harness should wait entirely — I do not think so, but it is a fair position and
this is the section to argue it in.
