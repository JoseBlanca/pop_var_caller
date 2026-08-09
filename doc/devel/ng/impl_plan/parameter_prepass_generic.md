# ng parameter pre-pass, the SNP/indel path (step 4) — implementation plan

**Status:** draft, 2026-08-06. The build order for the **generic half of step 4**: the
`parameter_estimation/` module, the vocabulary it adds to `types.rs`, the two keyed
accumulators, the fitting machinery both paths share, and the four numbers a sample emits —
a per-read-group error rate, its heterozygosity, its homozygous-non-reference rate, and its
inbreeding coefficient. Design is settled in
[`parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md) (spec) and
[`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) (types &
interfaces), on the shared framing of
[`parameter_prepass.md`](../spec/parameter_prepass.md) and under the shared arch docs
([step interfaces](../arch/ng_step_interfaces.md), [module layout](../arch/module_layout.md)).
This turns that design into build order; it is **not** a place for new design.

**What this step can be checked against, and why the answer is unlike the earlier plans'.**
Every ng plan so far has had one of two kinds of oracle, and this one has neither.

- **Parity with production**, where ng recomputes something production already computes: read
  filtering is drop-parity, read preparation is byte-parity with `process_read`, the STR
  generator is byte parity with production's tract tally. **Not available here, and not
  because production lacks an estimator — because it has a biased one.** Production returns
  `None` for a column carrying no alternative allele and its caller skips it
  ([`het.rs:146-148`](../../../../src/sample_summary/het.rs)), so the sites that are the
  majority of the genome and the strongest evidence there is about the error rate never enter
  the tally at all; and it classifies each remaining genotype before counting it. Those are
  the two biases this step exists to remove
  ([`parameter_prepass.md`](../spec/parameter_prepass.md) §2.1). **Agreeing with it would be
  the bug.**
- **A fixture whose right answer you can write down.** Where ng has no production oracle it has
  so far had this instead: the pileup locus generator's own plan says in as many words that
  beyond its parity class "there is no oracle, so the new behaviour is pinned by fixtures
  written to fail the *wrong* implementation"
  ([`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md)). That works
  because a locus's right answer is *legible*: lay out reads over a stretch, and what the loci
  and observations must be is something a person can state and a test can assert. **Not
  available here either.** Every number this step emits is the argmax of a sum over hundreds of
  millions of sites. There is no arrangement of reads whose fitted error rate a reader can
  write down, and no consumer that misbehaves visibly when it comes out wrong — which is the
  whole reason step 4 got a research note before it got a plan.

**So the oracle is a third kind: the estimator's bias, computed exactly.** Replace each cell's
observed count with the cell's probability under a known truth, maximise, and the answer is
what an infinite genome returns — a fixed number with no sampling noise in it, so "unbiased"
is decided rather than estimated. That is what the two harnesses in `examples/` do, and their
answers are in
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md).
Every milestone below is proven against them or against an identity, never against itself.

---

## Scope

**In:** `src/ng/parameter_estimation/` — `mod.rs`, `fitting/{mod.rs, mixture_weights.rs, ladder_scan.rs}`,
`generic/{mod.rs, accumulators.rs, coupled_fit.rs, depth_and_alt_reads.rs, depth_bins.rs,
fallback.rs, histogram.rs, noise_model.rs, read_group_error_rate.rs, runs.rs}`; the four constrained
newtypes step 4 adds to `types.rs`; `DepthBinEdges` and the cell table; the read-group and
windowed accumulators with their merge; `fit_mixture_weights`, the `NoiseModel` seam and the
profile scan; the coupled error-rate/frequency loop; the runs model; the fallback ladder and
`ParameterEstimationError`; both entry points.

**Out (each handed to a named later plan):**

- **The STR stutter histogram and its noise model** — spec
  [`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) is settled and its
  architecture is **in draft**, so its plan waits on that settling. It is the **second
  implementor of `NoiseModel` and the second consumer of `fitting/`** (arch §4), which this
  plan builds and it reuses — which is the thing that will show whether the seam was cut in
  the right place. Two changes it asks of `fitting/` are already known and neither is a
  rewrite: `fit_mixture_weights` widens past three genotypes, and a multi-start maximiser
  lands beside `fit_by_profile_scan` rather than replacing it.
- **The two censuses** — spec
  [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md); **no
  architecture document yet**.
- **The cohort gather** — spec
  [`parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md); **no architecture
  document yet**. `CohortEstimator` and the `SampleSummary` assembly
  ([`ng_step_interfaces.md`](../arch/ng_step_interfaces.md):351) belong to it.
- **How many samples are walked at once** — the multiplier that decides peak memory, and
  unchosen ([`parameter_prepass.md`](../spec/parameter_prepass.md) §6, spec §9). It is a
  property of the driver, not of these accumulators; it belongs to the plan that writes the
  per-sample walk.
- **Where `PloidyMap` gets its answer** — a flag, a BED or a per-contig default is not this
  module's business (arch §3). This plan builds the trait and a constant implementation,
  which is today's behaviour exactly.
- **The reference-bias term in place of the `½`** — spec §11.3 is open, and §8 now scopes what
  it would buy: heterozygosity on shallow cohorts and nothing above 20 reads a site (research
  note §5). Adopting it later lengthens the scan's ladder and changes no signature, because
  `NoiseModel::NoiseParams` is already an associated type carrying three parameters on the STR
  path (arch §4.2).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The scoring rule and the fits (Milestones D
  and E) are built and proven against the harnesses' exact answers **before** a single locus
  is read. What must be right is the mathematics; the accumulator only has to hand it the
  right counts.
- **Simplest source of data first, as the oracle for the next.** A histogram filled directly,
  cell by cell, is the fixture for every fit (E, F2). Only after the fits are proven does the
  locus stream feed one (F3). An accumulator bug and a fit bug then cannot hide each other.
- **Verify against ground truth that is not this code.** The harnesses' bias is computed
  exactly rather than simulated, so "matches the harness" is a real assertion and not a
  tautology. Where an identity is available it is asserted directly instead (arch §9, spec
  §12.8).
- **Isolate the steps whose failure is silent, and say so.** Most of this module fails loudly
  — a panic, a test. **Six do not**: the depth ladder (A4), the depth a cell is scored at (B3),
  the depth cap's draw (C2), the multi-library scoring rule (D2), the coupled loop (E2) and the
  runs model (E3). Each of those returns a plausible number nobody can check, and each has
  already produced a wrong one during the measurement work. They land as **their own commit
  with their oracle green before and after**, so `git bisect` can find one if a parameter later
  moves. They are marked **own commit, do not bundle** below; no other step is.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); a native
  host build at completion.

## Preconditions (already in place)

- **The generic locus stream exists and is complete.**
  [`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md) is finished
  (all seventeen steps), so `SampleLocusObservationsIterator`
  ([`locus_generation/mod.rs:701`](../../../../src/ng/locus_generation/mod.rs)) yields
  `LocusKind::Generic` loci with `complete_observations()` and `region` — the only two things
  this step reads from a locus (arch §3, contract).
- **`types.rs` carries the shared vocabulary** it builds on: `ContigId`, `GenomeRegion`, `Bp`,
  `ReadGroupId`, `LogProb`, `MismatchFraction` (the checked-newtype shape to copy) and
  `DomainError`, whose doc already names `InbreedingF` as expected
  ([`types.rs:268`](../../../../src/ng/types.rs)).
- **Both research harnesses run and their numbers are recorded.**
  [`examples/ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs)
  and [`examples/ng_inbreeding_harness.rs`](../../../../examples/ng_inbreeding_harness.rs),
  written up in the research note. They are the oracle for Milestones D and E and must be
  green before either starts.
- **No production dependency.** `src/ssr/` and `src/pileup/` are frozen;
  `sample_summary/het.rs` and `var_calling/posterior_engine.rs` are **shape to copy, not code
  to call** (arch §7).

---

## The steps

### Milestone A — vocabulary and the local types (types, no logic)

**A1. Scaffold the `parameter_estimation/` module tree.**  ✅
`mod.rs`, `fitting/{mod.rs, mixture_weights.rs}`, `generic/{mod.rs, depth_and_alt_reads.rs,
histogram.rs, runs.rs}`, each with its `#[cfg(test)]` block; wire `pub mod
parameter_estimation;` into `ng/mod.rs`. The folder split is the project rule that the shaping
of data and the mathematics on it never share a file. A fifth file under `generic/`,
`depth_bins.rs`, arrived with A4 and is now in the architecture's module table — the binning
rule fits neither side of that split, and all three A4 reviewers judged its own file right.
*Source:* arch §Module home, [module layout](../arch/module_layout.md).

**A2. Extend `types.rs` with the four constrained newtypes.**  ✅
`ErrorRate`, `GenotypeFrequency`, `InbreedingF`, `Ploidy` — each with a private field, a
`try_new` returning `DomainError`, and `.get()`, copying `MismatchFraction`'s shape
([`types.rs:243`](../../../../src/ng/types.rs)). Four types and not one shared `Probability`:
they are all fractions in `[0, 1]`, so one type would let an inbreeding coefficient be handed
to something expecting an error rate and compile. `Ploidy` rejects zero, because the
likelihood divides by it. Unit tests: boundary values accepted, out-of-range rejected, `Ploidy
= 0` rejected. *Source:* arch §2.1.

**A3. The step-4-local scalars and the error-rate ladder.**  ✅
`WindowIndex`, `INBREEDING_WINDOW_BP = Bp(100_000)`, the three `ERROR_RATE_LADDER_*_PHRED`
constants and `error_rate_ladder()`. Unit test: 161 rungs, ascending, spanning Phred 10 to 50,
adjacent rungs a ratio of `10^0.025` apart. No `Phred` newtype — a second log-scaled
probability type beside `LogProb` would make a base mix-up a plausible wrong number instead of
a compile error. *Depends:* A2. *Source:* arch §2.1.

**A4. `DepthBin` and `DepthBinEdges` — the ladder itself.**  ✅ **Own commit, do not bundle.**
Exact integers to 8, then eleven geometrically widening bins to a cap of 124 — twenty bins,
583 cells. `bin_for`, `row_start`, `depth_range` (a `RangeInclusive`, so an off-by-one cannot
silently mis-size a row), `cell_count`, `bin_count`, `max_depth`. **The silent failure this
isolates:** the edges are a correctness parameter, not a memory one — sixteen bins at the same
cap biases the error rate by 0.55 rungs and the homozygous-non-reference rate by 1.8%, against
0.05 rungs and 0.3% for twenty, and nothing downstream would show it (research note §4.3).
*Oracle:* assert the bin tops are exactly `10, 13, 17, 22, 28, 36, 46, 59, 75, 97, 124` above
the exact region, that `cell_count()` is 583, and that `bin_for` is monotone and total over
`0..=124`. *Depends:* A1. *Source:* arch §2.2, spec §4, research note §4.3.

**A5. The output types.**  ✅
`Provenance`, `Estimate<T>`, `SampleRates` with its two diploid accessors,
`GenericSampleParameters`, `FitTermination`, `ScanResult<P>`, `CoupledFit`, `RunsModelStarts`
(defaults `separations = [0.05, 1/3, 0.75]`, `implied_f = [0.05, 0.5, 0.75]`), `RunsModelFit`,
`StartOutcome`. Types only. Unit test: `SampleRates::observed_heterozygosity()` is `None`
above ploidy 2 and the frequencies sum to one. *Depends:* A2. *Source:* arch §2.4, §5.2, §5.3.

**A6. `ParameterEstimationError` and the fit floors.**  ✅
The four variants — `GenotypeFrequenciesNotFittable`, `InbreedingNotFittable`,
`InbreedingStatesNotSeparated`, `Domain` — with `MIN_SITES_TO_FIT = 10_000`,
`MIN_WINDOWS_TO_FIT_INBREEDING = 3_000`, `DEFAULT_ERROR_RATE = 0.001` and
`MAX_COUPLED_FIT_ITERATIONS = 20`. Test: each message names the sample and the number that was
too small. `InbreedingStatesNotSeparated` exists because a failed search and an outcrossing
genome leave identical fitted values, so returning zero there is the one way this estimator
produces a confident wrong number. *Depends:* A2. *Source:* arch §5.4.

> **Checkpoint A:** the vocabulary compiles; the ladder's edges and cell count are pinned by
> test against the measured ladder; every constrained newtype rejects what it must. Pause for
> review.

### Milestone B — the cell table (storage, no loci)

**B1. `SiteKey`, `DepthAndAltReads` and `CellCounter`.**  ✅
Two arms and no third: `Attributed { depth_bin, alt_by_group }` for at most
`MAX_ATTRIBUTED_ALT_READS = 4` alternative reads, `Pooled { depth_bin, alt_reads }` above it.
`alt_by_group` in read-group order so the key is canonical. `CellCounter` implemented for
`u32` and `u64` only — the trait exists to make the widening at the fold explicit in the type.
Unit test: two sites differing only in the order their read groups are listed produce the same
key. *Depends:* A4. *Source:* arch §2.2.

**B2. `DepthAltHistogram<C>` — storage, `add_site`, `cells`, the two counters.**  ✅
The flat ragged `counts` and `depth_sums` located through `edges.row_start`, the sparse `fine`
map for the attributed arm, the `Arc<DepthBinEdges>` handle. `add_site` derives the bin from
the exact depth it is handed and adds that depth to the cell's running sum. `cells(ploidy)`
materialises once — a `Vec` and not an iterator, because the profile scan re-walks it 161
times. `total_loci` counts entries; `total_covered_positions` counts reference bases, and the
two differ because a generic locus can be widened to an indel's reference span. Unit tests:
a hand-built table's rows are the right widths; `cells()` is stable in order across runs.
*Depends:* B1. *Source:* arch §2.2.

**B3. `mean_depth_in_cell`.**  ✅ **Own commit, do not bundle.**
The mean of the exact depths that landed in **this cell**, from its own depth sum. **The
silent failure this isolates:** taking the mean over the whole *bin* instead charges 0.3% of
sites a negative number of reference reads, and the fit then lands 5.2 rungs below the true
error rate and 29% below the true homozygous-non-reference rate — bounded, so
`argmax_at_ladder_end` never fires, and there is nothing on the outside to see (research note
§4.5). *Oracle:* the identity `cell.alt_reads() ≤ mean_depth_in_cell(cell)` asserted at every
cell of a table built to hold depths 100–124 with alternative counts up to 124 — which is the
exact case the per-bin mean fails — plus a per-bin-mean unit test showing that same table
violating it, so the assertion is proven to bite. *Depends:* B2. *Source:* arch §2.2, spec
§12.10, research note §4.5.

**B4. `merge` and `fold_windows_of_one_ploidy`.**  ✅
Element-wise integer addition on the pooled table and a key-wise sum on the attributed map,
panicking unless the two histograms hold the same edges object (`Arc::ptr_eq` — a proof, not a
length comparison). `fold_windows_of_one_ploidy` — named `whole_sample_histogram` when this
plan was written, renamed on the owner's call (2026-08-06) because the ploidy restriction
cannot live in the signature and "whole sample" reads as *all* of it — folds the windows for
one ploidy and **widens both counters to `u64` here and only here**. The depth sum is the one that forces it: folded
over a human genome the site count reaches 3.1 × 10⁹ against a `u32` ceiling of 4.29 × 10⁹ —
close, but inside — while the depth sum reaches 3.1 × 10¹¹, **seventy-two times over**. A fold
that widened the site counts and left the depth sums alone would wrap the very quantity
`mean_depth_in_cell` exists to hold, which is B3's failure by another route. Unit test: a table split arbitrarily in two and merged equals the
unsplit one, cell for cell, in either merge order. *Depends:* B3. *Source:* arch §2.2, §3.

> **Checkpoint B:** the table stores, merges and reports; the depth a cell is scored at is
> proven never to fall below its own alternative count. Pause for review.

### Milestone C — one locus → one cell (data shaping)

**C1. `count_whole_site` and `count_by_read_group`.**  ✅
The only place that decides what counts as an alternative read — which is why it is its own
file and not a method on the locus type. **Complete witnesses only**
([`locus_generation/mod.rs:134`](../../../../src/ng/locus_generation/mod.rs)): a read that
spanned part of the locus witnessed neither allele at the positions it missed.
`reads_without_observation` does **not** enter the depth — those reads showed nothing, and
counting them would assert they showed the reference. Unit tests over hand-built
`SampleLocusObservations`: a one-base locus, a locus with a partial witness, a locus where
every read is a non-witness. *Depends:* B1. *Source:* arch §2.3.

**C2. The depth cap — subsample, do not rescale.**  ✅ **Own commit, do not bundle.**
A site deeper than `max_site_depth(edges)` keeps 124 of its reads and counts the alternative
ones among them, seeded from the locus position so a region-sharded walk and a single-threaded
one keep the same reads. Fires in `count_*`, before the pair is built, so the depth recorded
and the depth the alternative count belongs to are the same number. **The silent failure this
isolates:** the two tempting shortcuts are both quietly wrong — rescaling and rounding to
nearest reverses the bias's sign at the depth where a lone alternative read stops surviving
the round, and a stochastic round fixes the mean while making the spread four times too
narrow, which the fit reads as a cleaner read group than the data supports. *Oracle:* over
many seeds the kept alternative count is hypergeometric — mean and variance match the closed
form — and the same locus position gives the same draw on every run. *Depends:* C1, A4.
*Source:* arch §2.2 (`max_site_depth`), §2.3.

**C3. `GenericAccumulators`, `add_locus`, `AccumulationCounts`.**  ✅
The two keyed collections (`BTreeMap` and not `HashMap`: the runs model reads windows in
genome order, and every fit is a floating-point sum over cells, which is not associative);
`InbreedingMode`, which drops the window key when `F` is supplied and collapses the object
from 37 MB to a few kB; `PloidyMap` with a constant implementation. `add_locus` borrows and
passes the locus on untouched, ignores a non-generic `kind`, and tallies rather than repairs:
`loci_with_upstream_subsample`, `reads_without_observation`, `sites_subsampled_to_cap`, and
`loci_overlapping_previous`, **which must read zero**. The overlap check keeps the *span each
shard covered* as a list that `merge` concatenates and sorts once at the end — collapsing it
to a first-start/last-end pair per contig reports a false overlap when three contiguous shards
merge out of order, and a counter that must be zero may not have false positives. Unit tests:
a non-generic locus changes nothing; two overlapping loci are counted; three shards merged in
every order give the same counters. *Depends:* C2, B4. *Source:* arch §3.

> **Checkpoint C:** loci reduce to cells, the cap is exact, and sharded accumulation is proven
> order-independent including the counter that must be zero. Pause for review.

### Milestone D — the fitting machinery (the mathematics, no loci)

**D1. `fit_mixture_weights` — the concave climb.**  ✅
Given each cell's per-genotype likelihood and a weight per cell, the genotype frequencies that
best explain the table. Expectation-maximization is a reasonable default and nothing depends
on it being EM. **Convergence failure is a bug, not a data condition** — the surface is
concave, so it is asserted in tests rather than propagated as a flag no consumer would read.
Not used by the runs model, whose two states are a constrained parameterisation rather than a
free point on the simplex. Unit test: recovers known frequencies from a hand-built table from
any interior start. *Depends:* A5. *Source:* arch §4.1, spec parameter_prepass §3.1.

**D2. The generic `NoiseModel` — §5.1's closed form.**  ✅ **Own commit, do not bundle.**
`ln L(cell | θ)` summing over the split the key forgot rather than inventing a per-library
depth: one multinomial over `G + 1` categories, one "alternative from library g" per library
and one pooled "showed the reference". **The silent failure this isolates:** the plug-in an
earlier draft used — give each library `n̂_g = w_g·n` — is not a probability over the cell
space and reports heterozygosity 68% high at three reads on two libraries with the *same*
error rate, and it does not shrink as data accumulate. *Oracle:* the three algebraic identities
of spec §12.8, as unit tests, none needing a fit — the rule sums to one over the cell space at
any parameter values; no cell is charged a negative count of reference reads; and with every
library's error rate equal it reproduces the exact per-library likelihood to floating point.
Then a fourth: on the cell space of one of the harness's worlds, this implementation's
`ln L` agrees with
[`ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs)'s
`ln_component_attributed` to floating point. **Any later change to this expression re-runs all
four first.** *Depends:* B3. *Source:* arch §5.1, spec §1, §12.8.

**D3. `fit_by_profile_scan` and the rail flag.**  ✅
Step through the ladder, climb to the best frequencies at each rung, keep the best-scoring
rung. **Every rung is scored — no early exit**, because nobody has shown the curve has a
single hump. Takes a slice of cells and not a histogram, which is what makes the "shared with
the STR path" claim true rather than aspirational. Ploidy travels with each cell, because one
error rate is fitted per read group across every ploidy that group covered while the
frequencies are climbed once per ploidy. `ScanResult` carries `argmax_at_ladder_end`, and it
is the only thing between a railed fit and a plausible-looking number. Ties resolve to the
lower error rate, stated so two implementations cannot differ. Unit tests: a table generated at
a known rate recovers that rung; a table generated outside Phred 10–50 sets the rail flag.
*Depends:* D2, D1. *Source:* arch §4.2.

> **Checkpoint D:** the scoring rule passes all four identity checks and the scan recovers a
> known rate from a synthetic table. Nothing has read a locus yet. Pause for review.

### Milestone E — the four fits

**E1. The per-read-group error rate.**  ✅
A scan over `ReadGroupHistograms`, once per read group, keeping only `ε`. **Not
`fit_by_profile_scan`** — that climbs its own frequencies at every rung, which is a different
estimator and the one never measured. E1 is the `ε` half of E2's alternation, so it scores
every rung at the genotype frequencies it is **handed**, one shared set across the read
groups, and climbs nothing (owner's call, 2026-08-07; the harness's own
`fit_eps_on_read_group(space, freqs)`). That is a **sibling** of the profile scan rather than
a mode of it, because a scan at fixed frequencies is not a profile likelihood.
*Depends:* D3, C3. *Source:* arch §5.1, spec §3, §5.1.

**E2. The coupled loop.**  ✅ **Own commit, do not bundle.**
Alternate: the frequencies from the whole-sample table at the previous rates, then each read
group's rate from its own table at **those** frequencies and without re-climbing them —
capped at 20 iterations, keeping the **best-scoring** iterate and reporting `FitTermination`.
(An earlier draft had the two blocks the other way round and the rate step re-climbing; the
order is a phase and changes no fixed point, but the re-climbing is a different estimator.) **Stop when every read group's winning
rung is unchanged** — the scan returns a rung index, so "moves by less than one rung" and "does
not move" are the same condition and only the second is testable. **The silent failure this
isolates:** this is a fixed point of two estimating equations rather than a climb on one
objective, so a wrong alternation converges to a plausible wrong pair and reports success; a
loop oscillating between two adjacent rungs would satisfy a movement tolerance forever.
*Oracle:* from a start at three times the true rates and half the true frequencies, the fixed
point is the truth in the harness's 25 worlds — error rates to 0.000 rungs and both
frequencies to 0.000% (research note §2.6). And **at one read group the alternation must reach
the profile scan's answer**: with one library the two tables are the same table, so both
procedures converge to the same joint maximum — each block being an exact maximisation of one
objective. (It does **not** terminate after one iteration, which an earlier version of this
oracle asked for: that is true of a profile scan and false of coordinate ascent. What the
difference costs is iterations, not answers.) *Depends:* E1. *Source:* arch §5.2, spec §5.1.

**E3. The runs model — a two-state HMM over windows.**  ✅ **Own commit, do not bundle.**
Each state its own three genotype frequencies, fitted freely, with the ordering constraint
`h << Hout` applied by relabelling after the fit; both transition rates fitted per base; the
emission a **sum over the window's cells**, never a per-window heterozygote count; the chain
covering every window of a contig, absent ones included as empty, restarting at each contig
boundary. `F` is the coverage-weighted posterior occupancy — weighted by
`total_covered_positions()`, not by loci — and not the transition rates' ratio; the two differ
by 3.5% to 11% on a finite genome. **The silent failure this isolates:** starts that disagree
only about `F` are not a spread. Nine starts spanning the state separation return `F` = 0.2634
where five sharing one separation guess return `F` = 0.0000, converged and silent, on the same
genome (research note §3.4). Climb from every start in `RunsModelStarts`, keep the
best-scoring, **report them all in `starts_tried`**, and return
`InbreedingStatesNotSeparated` — never zero — when no start left mass on both states.
*Oracle:* [`ng_inbreeding_harness.rs`](../../../../examples/ng_inbreeding_harness.rs) — `F`
recovers a drawn genome's **realised** autozygous fraction to four decimal places at true
values 0.05, 0.15, 0.30 and 0.60, and a floor of false heterozygotes at up to five times the
real rate does not move it. Score against the realised fraction and never the nominal one: a
finite genome does not have the `F` its rates imply, and comparing against the nominal value
reads sampling as bias. **It takes the per-read-group error rates rather than one pooled
rate**, so it follows E1 — a site with few alternative reads keeps which library each came
from, and each must be weighed against its own library's rate. **It does not use
`fit_mixture_weights`**: its two states are a constrained parameterisation over `(f, h)`, a
surface inside the simplex rather than a free point on it, and the concavity that makes D1
safe does not transfer to a curve. *Depends:* E1, C3, B3. *Source:* arch §5.3, §4.1, spec
§6.1, §6.5.

**E4. The fallback ladder and the floors.**  ✅
Fitted here → borrowed from the sample's other groups → supplied → defaulted, each carrying
its `Provenance`. `MIN_SITES_TO_FIT` gates the first; below
`MIN_WINDOWS_TO_FIT_INBREEDING` the runs model **fails rather than emits**, because inbreeding
is the parameter that differs most between an outcrosser and a selfing landrace and the
cohort's diversity divides by `1 − F`. `RunsModelFit::resolution` reports the noise floor at
this run's window count — about 0.01 at tomato's 8,004 windows — so a consumer can tell *no
runs detected* from *a small autozygous fraction*. Unit tests: a thin read group is marked
`Borrowed`; every group thin gives `Defaulted`; 1,000 windows returns
`InbreedingNotFittable`. *Depends:* E3. *Source:* arch §5.4, spec §6.1.

> **Checkpoint E:** all four parameters are fitted, and each is proven against the harness
> answer it has to reproduce, not against itself. Pause for review.

### Milestone F — the entry points and end to end

**F1. The two ways in.**  ☐
`GenericEstimationConfig`, `estimate_generic_parameters(loci, config)` for a caller with
nothing else to do with the stream, and `GenericAccumulators::estimate(config)` for one that
drove the accumulator itself. The first is the second over an accumulator fed by the stream,
so the two cannot diverge. A `LocusGenerationError` in the stream is fatal and propagates: the
loci a walk failed to produce are missing evidence, not zero evidence, and a rate fitted over a
truncated genome is wrong in a way nothing announces. *Depends:* E4. *Source:* arch §1.1.

**F2. Recovery from a directly-filled accumulator — no reads, no reference.**  ☐
Fill a histogram cell by cell from known parameters and refit. At ploidy 2 **and 4**, and
**at 3 reads a site and at 300×**: at tomato's depth every site sits in a one-per-depth bin,
so a binning fault is invisible below about 100×, and the deep arm is the only one that
exercises the cap of C2. **The tolerance is one rung of the error-rate ladder**, which is the
resolution the design already argues is finer than a caller can feel
([`parameter_prepass.md`](../spec/parameter_prepass.md) §3) — not a number chosen here. For
the two frequencies the measured binning bias of the adopted ladder is 0.3% (research note
§4.3), so assert 1% relative: loose enough that binning alone cannot fail it, tight enough
that a real fault cannot pass. *Depends:* F1. *Source:* arch §9, spec parameter_prepass §10.1.

**F3. The identities that need no simulated truth, on both real cohorts.**  ☐
Three assertions on the tomato CRAMs and the HG002 alignments as they stand, all of which hold
by construction and none of which needs a truth set: the read-group histogram equals the
windowed one folded over its windows, cell for cell, on a single-library sample (which is every
sample in both cohorts); one sample walked in one region and in many gives identical
histograms; and `adjustments().loci_overlapping_previous` is zero. The last is a bug report
against locus generation rather than something this unit absorbs. *Depends:* F2. *Source:*
arch §9, spec §12.6.

> **Checkpoint F:** step 4's generic path runs end to end on real alignments, and the three
> structural identities hold on both cohorts. Pause for review.

### Milestone G — the anchors against real data

**G1. Model-free values from the GIAB truth set.**  ☐
The error rate as non-reference bases over total bases at truth homozygous-reference
positions — a count, no model and no fit; heterozygosity as the truth het count over the
confident regions' length; the homozygous-non-reference rate as the 1/1 count over the same.
The truth set is `benchmarks/giab/all_bench_regions/`, the whole-genome HG002 v4.2.1 small
variant VCF and its confident BED — **not** `benchmarks/ssr_hg002/`, which is the tandem-repeat
benchmark and routes to the STR path. These **bound** the fitted values rather than pin them,
because confident regions are the easy regions: a fitted error rate below the model-free one
on easy regions is an unambiguous bug. *Depends:* F3. *Source:* arch §9.

**G2. Invariance to coverage, which needs no truth at all.**  ☐
Fit all four parameters on one HG002 alignment downsampled to 300×, 30×, 10× and 3×,
**restricted to the confident BED at every rung** or the arms compare different site sets
rather than different depths. Same genome, so an error rate (per read) and three properties of
a genome must all come out flat; any slope is bias and its sign names the mechanism. The
tomato equivalent is free: plot each sample's fitted heterozygosity against its mean depth
across the cohort, where a biological quantity has no business correlating with library yield.
*Depends:* G1. *Source:* arch §9.

**G3. `F ≈ 0` on HG002 — blocked on data, and the blockage is named here rather than
discovered later.**  ☐
HG002 is not consanguineous, so `F` must come back at the noise floor. **This is the one
anchor a region subset cannot carry.** The alignments in `benchmarks/giab/per_sample/bam/` are
100 randomly selected regions, about 1,200 windows, and at that count a genome generated with
**no runs at all** returned `F` averaging 0.23 and 0.84 on one seed of eight (research note
§3.6) — so restricting the anchor there does not weaken it, it voids it. The other three
parameters are fine on a subset because they are per-read or per-site rates. **This step needs
a whole-genome HG002 alignment fetched, and it has no substitute.** It is a data question, not
a design one; until it is resolved this checkbox stays open and `F` ships with no anchor
against real data. *Depends:* G2. *Source:* arch §9.

> **Checkpoint G:** three of the four parameters are anchored against something that is not
> our own model. `F` is anchored only if the whole-genome alignment was obtained; if it was
> not, that is recorded here as the one gap step 4 ships with. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the ladder's bin tops and 583-cell count pinned by test against the measured ladder (research note §4.3); every constrained newtype rejects out-of-range and `Ploidy = 0` |
| B | `alt_reads ≤ mean_depth_in_cell` asserted at every cell of a table spanning depths 100–124, **and** the per-bin mean shown to violate it on the same table; merge exact and order-independent |
| C | the cap's kept alternative count is hypergeometric in mean and variance over many seeds, and reproducible from the locus position; three shards merged in every order give identical counters, `loci_overlapping_previous` included |
| D | the four identities on the scoring rule — sums to one, no negative reference count, exact at equal error rates, and agreement with `ng_multilib_key_harness.rs` to floating point; the scan recovers a known rung and flags a railed one |
| E | the harnesses' own answers: the coupled fixed point is the truth in 25 worlds from a deliberately wrong start; `F` recovers a drawn genome's realised autozygous fraction to four decimal places and does not move under a false-heterozygote floor of five times the real rate; a start set sharing one separation guess is proven to return `InbreedingStatesNotSeparated` rather than zero |
| F | recovery from a directly-filled accumulator at ploidy 2 and 4 and at 3 reads and 300×; the two histograms equal cell for cell, sharded equals single, and the overlap counter is zero — on the tomato and HG002 cohorts as they stand |
| G | model-free counts from the GIAB truth set bound the fitted values; the four parameters are flat across a 100-fold coverage sweep on one genome. **`F ≈ 0` is blocked on a whole-genome HG002 alignment and has no substitute** |

## Out of scope (next plans)

- **The STR stutter path of step 4** — spec
  [`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) is settled and its
  architecture is in draft; its plan follows that. It is `fitting/`'s second implementor, and
  the two changes it needs there are additive: a `fit_mixture_weights` that widens past three
  genotypes, and a multi-start maximiser beside the profile scan.
- **The two censuses** — spec
  [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md); needs an
  architecture document.
- **The cohort gather, `SampleSummary` assembly, and `CohortEstimator`** — spec
  [`parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md); needs an architecture
  document. It owns the diversity, the frequency spectrum, and the contamination split.
- **The per-sample walk and its concurrency** — including how many samples are in flight, the
  multiplier on 37 MB per tomato sample and 145 MB per human one that decides peak memory
  (spec §9). A driver decision, not an accumulator one.
- **The memory measurement** that replaces spec §9's arithmetic
  ([`parameter_prepass.md`](../spec/parameter_prepass.md) §10.6), now pricing 583 cells a
  window rather than the 465 an earlier draft assumed.
- **What a site deeper than the cap costs** — the subsampling rule of C2 is specified and
  implemented here, but no harness measures its bias, because no world in either harness
  reaches a cap (research note §4.6). It fires on samples above ~124× — HG002, never tomato —
  and it is the one depth mechanism shipping without a measurement. Arch §8 carries it as the
  single remaining `OPEN:`.
