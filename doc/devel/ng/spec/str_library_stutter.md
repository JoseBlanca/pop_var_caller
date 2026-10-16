# Each library's own stutter: a multiplier fitted jointly with the repeat-class rates

**Status:** design, 2026-10-10, for the owner's review. Decided with the owner on 2026-10-10 (§1.3).
Build order: [`../../implementation_plans/str_library_stutter.md`](../../implementation_plans/str_library_stutter.md).
It amends the repeat-tract fit ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.2), the
curves across repeat count only in what they are fitted to
([`str_slippage_level_curve.md`](str_slippage_level_curve.md)), calling's stutter lookup
([`read_likelihoods.md`](read_likelihoods.md) §4.4) and the parameters file
([`parameters_file.md`](parameters_file.md), format version 3).

> **The fit this design sits on was replaced, 2026-10-16; the design waits to be rebuilt on its
> successor.** §3 is written around the climb: each repeat class fitted on its own by an alternating
> climb, restarted from its last answer every round (§3.1); a sample's genotype taken from that
> climb's likelihood, `ln_tract` (§3.2); a class read from a subset of samples (§3.5); and the curves
> drawn through the classes' answers, then blended with them (§3.6). Penalised expectation
> maximisation (EM), [`str_penalised_fit.md`](str_penalised_fit.md), replaced that fit; it fits
> every class of one motif length together with each number a smoothed curve along the repeat count,
> and has no subsets and no blend. Commit `3198d55d` deleted the climb, `ln_tract`, the subsets, the
> blend and this design's prototype (`ssr_fit/library_stutter.rs`). The new spec defers library
> multipliers to this document, to be rebuilt on the new fit (its §1.4 and §6).
>
> **Still true:** the problem and its measurements (§1.1), the owner's three decisions (§1.3), the
> model — a library's level is its class rate times its multiplier, scaled so the multipliers'
> weighted geometric mean is one (§2) — the pull toward one (§3.3), and how calling uses a
> multiplier (§4). **Still in the code:** each read group kept apart as a library in the evidence, at
> a multiplier of one (`gather_strata`), and `library_level` with its ceiling. **No longer right:**
> §3.1, §3.2's source of the genotype, §3.5 and §3.6, which name code that is gone; and §5's format
> version, since version 3 of the parameters file was taken by the new fit's companion,
> [`str_fit_to_calling.md`](str_fit_to_calling.md), without this table.

---

## 1. What this is

### 1.1 The problem

**The repeat-tract fit gives every library the same stutter.** `estimate-parameters` puts every
read group into one slippage group (`every_read_group_pooled`, `src/run/census_fit.rs:288`), so it
fits one slippage level for each repeat class — a motif length and a reference repeat count, the
fit's *stratum* — and calling scores every library's reads with it. The pooling was recorded as a
choice ("production pools every read group into one slippage group", [`fit_precision.md`](fit_precision.md)
§4.1) and never revisited.

**Libraries do not stutter alike.** On 51 tomato accessions the share of stuttered reads ran from
0.37% to 1.98% between samples, about fivefold, far beyond what the read counts allow by chance, and
partly following each library's duplicate rate, a sign of PCR
(`str_stutter_benchmark.md` §3.3, a spec on branch `tract-ceiling`). The owner
expects libraries' stutter levels to be very different in general.

**What pooling costs, measured on the simulated stuttering benchmark** (branch `tract-ceiling`,
report `str_stutter_benchmark_first_pass_2026-10-10.md`): each simulated library stutters at its
level's rate times its own multiplier, drawn log-normal with a standard deviation of 0.4 on the log
scale. In the 50-sample PCR cohort the multipliers run from 0.33 to 3.45, and 27 of the 50 libraries
stutter more than 25% away from the cohort's mean. The pooled fit recovers the mean well — its level
divided by the true mean rate has a median of 1.01 over 91 classes at 30× — but not any one library.
Calling with each library's true stutter called more cells right than calling with the fit in all 16
of the benchmark's runs, by 0.1 to 6.1 points.

### 1.2 What changes

**Each library gets one stutter multiplier, and its slippage level in a repeat class is the class's
rate times that multiplier.** The class rates and the multipliers are fitted together, so the class
rate is what the libraries share once each library's own level is accounted for — not the level of
whichever kind of library is most common in the cohort.

- **A library is a read group**, the grain the census already keeps distinct across the cohort (a
  cohort in which two samples claim the same read group is refused, `census.rs:2120`).
- **Only how often a read slips varies by library.** Which way it slips (the shorter share) and how
  far (the fall-off) stay one per repeat class. *Assumption, not measured*: the simulator varies only
  the rate, and the tomato spread was measured on the rate.
- **The substitution rate inside a tract stays pooled**: one rate a repeat class, counted over every
  library and copied to each library's row (`census_fit.rs:376`). Out of scope (§8).

### 1.3 Decided with the owner (2026-10-10)

1. **Not a shared stutter first and a multiplier after.** "Do not create a shared stutter, I don't
   want to risk biasing towards the most common stutter observed in the libraries." The class rates
   are refitted with the multipliers in place, alternating until the multipliers settle (§3).
2. **Thin evidence is pulled toward a multiplier of one**, by a few pseudo-reads, as the per-locus
   refit pulls a locus toward its class's rate (§3.3).
3. **One multiplier per library, not a full set of rates per library.** A library at 3 reads a
   position has too few slipped reads to fit its own rate in each of about 90 repeat classes; one
   number per library is well determined.

### 1.4 Vocabulary

- **Repeat class** (the code's *stratum*): a motif length and a reference repeat count.
- **Class rate** `Lₛ`: the slippage level of repeat class `s` for a library whose multiplier is one.
- **Multiplier** `mₗ`: library `l`'s stutter relative to the class rate. A read of library `l` over a
  tract of class `s` slips with probability `min(mₗ · Lₛ, LEVEL_CEILING)`.
- **Slippage group**: a set of libraries the run declares share one set of class rates. Today, and in
  this design, there is one. The multipliers live inside a slippage group (§2.2).

---

## 2. The model

### 2.1 A read's slip

A read crossing a tract of class `s`, from library `l`, reports a length different from its allele's
with probability `λₗₛ = min(mₗ · Lₛ, LEVEL_CEILING)`; given that it slipped, the direction and the
distance follow the class's shorter share and fall-off, as today
(`Slippage::read_probabilities`, `ssr_fit.rs:111`). `LEVEL_CEILING` = 0.95: a level of one would say
no read reports its allele, and the kernel's same-length share is floored at 0.01 downstream
(`StutterModel::sanitized`, `alignment/stutter.rs:263`).

### 2.2 What fixes the scale

`mₗ · Lₛ` is unchanged when every multiplier in a slippage group is doubled and every class rate
halved. **The multipliers of a slippage group are scaled so that their geometric mean, weighted by
each library's spanning reads, is one**, and the class rates absorb the factor. So a class rate is
the level of a library of typical stutter, and a multiplier says how far a library is from typical —
which is what a reader of the file needs, and what the pull toward one (§3.3) is a pull toward.

### 2.3 One library in a slippage group

**Its multiplier is one, written `defaulted`.** With one library there is nothing to compare it
against: its own level is the class rate. This is the one-sample case and the case of a cohort whose
samples were each sequenced once from one library pool declared as one read group.

---

## 3. The fit

### 3.1 The alternation

1. **Round 0 is today's fit**: every multiplier one, every class fitted as now. *Its numbers equal
   today's to the bit*: with every multiplier one, scoring a library's reads at `1 · Lₛ` is scoring
   the pooled reads at `Lₛ`.
2. **Then, each round:**
   1. **The multipliers, with the classes held** (§3.2), over every sample — every library gets one,
      whether or not its sample was in a class's subset (§3.5).
   2. **Scale** (§2.2).
   3. **The classes, with the multipliers held**: every class refitted, each read scored at its own
      library's level, each climb starting from the class's last answer.
3. **Stop** when every multiplier moved by less than `SETTLED_FRACTION` (0.1) of its own standard
   error over the round — the rule the rest of the fit stops by ([`fit_precision.md`](fit_precision.md)
   §2) — or after `MAX_MULTIPLIER_ROUNDS` (soft, 5), and say which.

### 3.2 One library's multiplier, with the classes held

**Each read is shared between "slipped" and "read exactly" by how likely each is**, given the
sample's genotype at the tract as the classes' fit sees it. The expected number of slipped reads
`Sₗₛ` and of exactly-read ones `Nₗₛ`, summed over library `l`'s tracts in class `s`, are what one
pass over the evidence gives, per library and class. With the classes held, library `l`'s
log-likelihood in its multiplier is then

`Σₛ [ Sₗₛ · ln(m · Lₛ) + Nₗₛ · ln(1 − m · Lₛ) ]`

— concave in `m`, one dimension, solved by Newton's method inside `(0, LEVEL_CEILING / maxₛ Lₛ)`.
This is the EM step for the multiplier; repeated rounds of §3.1 are what make it the joint maximum.

**The genotype a sample carries at a tract** is the posterior the class's own likelihood gives:
the sample's reads under each genotype, times the genotype's probability under the tract's
allele-frequency integral and the sample's homozygote excess (`ln_tract`, `ssr_fit.rs:3884`).

### 3.3 The pull toward one (owner decision 2)

**`P` pseudo-reads that slipped at exactly the expected rate are added to every library**: the
objective gains `P · (ln m − m)`, whose maximum on its own is `m = 1`. A library with many slipped
reads is barely moved by it; one with few is held near one. `P` is soft: starting value 20, the
per-locus refit's own (`DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS`), to be set from the benchmark at 5×.

### 3.4 Its standard error

The curvature of the objective of §3.2 at its maximum, inverted: about `m / √(slipped reads)` where
the cap does not bind. Written beside the multiplier (§5).

### 3.5 With a class fitted on a subset of samples

A class can be fitted on a subset of the cohort's samples ([`fit_precision.md`](fit_precision.md)
§4.4). **The multipliers are always estimated from every sample**; the class's climb then reads its
subset's libraries at their multipliers. The subset's top-up rule stays on slippage groups, not
libraries — every library a group of its own would force every sample into every subset.

### 3.6 The curves and the blend

Unchanged in shape: the curves across repeat count are fitted to the class rates `Lₛ`, per slippage
group and period, and each class's rate is blended with its curve as today. A class taken from its
curve gets its rate from the curve, and each library its multiplier on top.

---

## 4. Calling

**A cell's frozen level is its class's rate times its library's multiplier, capped at
`LEVEL_CEILING`**, set where calling gathers a locus's stutter
(`TractScoringFits::gather_for_locus`, `calling/inference/repeat_tract_parameters.rs:338`). The
per-locus refit then works on that level unchanged: its expected slips and its pull-back target are
the library's own (`slippage_refit.rs`).

- **`--defaults`, and a file without multipliers**: every multiplier one.
- **A read group the file has no multiplier for**: one, counted beside the read groups with no
  slippage group in the run report (`RepeatTractFitsUsed`).

---

## 5. The parameters file, version 3

**A new table, `repeat_tracts.stutter_multiplier_by_read_group`**: one row a read group,
`{ read_group, multiplier = { value, warrant, observations, standard_error } }`, the shape of
`base_quality_calibration.by_read_group` ([`parameters_file.md`](parameters_file.md) §3.3).

- `warrant`: `fitted_here`, or `defaulted` for a library alone in its slippage group (§2.3); a person
  may write `supplied`.
- `observations`: the library's expected slipped reads over all classes.
- `FORMAT_VERSION` 3. **A version-1 or -2 file reads, with every multiplier one**; a version-2 build
  refuses a version-3 file with its existing message.
- The file's explanatory text says what the multiplier is, that a class rate is the level of a
  library of typical stutter (§2.2), and that `mₗ · Lₛ` is capped.

---

## 6. What moves

- **Every fitted number of the repeat-tract half can move** once a cohort has more than one library,
  and the calls with it. **Round 0 is today's fit**, so a cohort of one library moves nothing; the
  four-accession oracle (four libraries) moves, and the move is measured and explained before its
  checksums are re-recorded.
- **Time.** Each round refits every class, from its last answer. On kimura the repeat-tract fit took
  10 h 41 min cold; a round is warm and cheaper, by an amount to be measured. The multipliers'
  pass is one scoring pass over the evidence a round.
- **Memory.** The evidence keeps each library's counts apart rather than summing them into the
  slippage group's. With one library a sample — every sample in both benchmark cohorts — the size
  does not change.

---

## 7. How we know it works

1. **Drawn classes with known multipliers** (unit tests): a cohort drawn at known class rates and
   known per-library multipliers is fitted; every multiplier lands within two standard errors of its
   truth in about 95 cases in 100, at 3 and 30 reads, 2 to 50 libraries.
2. **Round 0 equals today's fit to the bit** (unit test), and a one-library cohort equals it after any
   number of rounds.
3. **The simulated stuttering benchmark** (branch `tract-ceiling`): the recovered multipliers against
   the simulator's, library by library; and calling accuracy with the multipliers against the pooled
   fit and against the true stutter, at 50 samples and 5× and 30×, and at one sample. **The bar**: the
   multipliers recover the simulator's to within their errors, and the fit closes most of the gap to
   the true stutter at 50 samples.
4. **The oracle cohort and the cross-platform checksums**: moved, measured, explained.

---

## 8. Deferred

- **The substitution rate per library** — pooled today and copied to each library's row. Home: this
  document's successor, if the benchmark or a real cohort shows libraries differ in it.
- **Per-library direction and fall-off** — only if a measurement shows libraries differ in them.
- **Slippage groups other than one** — declared by the run (a flag the parameters file already has a
  table for, `slippage_group_by_read_group`); the multipliers live inside each group.

## 9. Open questions

1. **`P`, the pull toward one** (soft, 20) — settled by the benchmark at 5×.
2. **`MAX_MULTIPLIER_ROUNDS`** (soft, 5) — settled by how many rounds the benchmark and kimura take.
