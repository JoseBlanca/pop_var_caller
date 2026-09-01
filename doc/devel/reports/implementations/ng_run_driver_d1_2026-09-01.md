# ng direct mode, step D1 — calling, joined to the merge

**Date:** 2026-09-01. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step D1. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md)
§3.1, §5.1, §8; [`../../ng/spec/cohort_merge.md`](../../ng/spec/cohort_merge.md) §3.3, §6.1, §6.3.
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §3.1, §3.2,
§3.4, §5. **Modules:** `src/ng/run/callers.rs`, `src/ng/run/mod.rs`,
`src/ng/run/cohort_merge/{build,serial,observation_cache,timing}.rs`,
`examples/ng_open_cohort_descriptors.rs`.

---

## What landed

**`AlignedFilesVariantCaller::call_cohort` — reads in, called loci out.** Every sample's
alignment files are walked at the merge frontier, the cohort's loci are assembled over the
analysed ground, and each locus is genotyped **where it is built**, before the next one is
closed. The three-call chain arch §3.2 sketches — `select_generic` → `shape_generic_locus` →
`LocusGenotyper::call_locus` — had no caller anywhere in the tree before this step; composing it
is what D1 is.

**The two owner decisions of Checkpoint C landed here.**

- **A run's walk tallies survive the merge.** `ObservationCache::into_sources` hands the
  per-sample readers back, so each walker's region accounting and its SNP/indel generator's
  counters are copied out before the walkers are dropped, along with the assembly-check outcome
  the run computed at construction. What is still out of reach is the per-read-group read-filter
  tallies, and the obstacle is not a missing accessor: at each contig boundary the retiring
  cursor's read-group counts are dropped rather than accumulated, so by the end of a walk every
  contig but the last has already lost them. Spec §8 requires a run to sum them at the end;
  reaching them is a change to the generator.
- **A run can set its locus generator's settings.** `AlignmentInputs` gains a sixth field,
  `locus_generator_settings`, checked at `open` with the other five refusals — before a file is
  opened, so a thousand-sample cohort is not opened to be told at its first locus.
  `RunError::LocusGeneratorSettings` was unreachable and now fires; its message is the cause's
  alone, since a wrapper would put a sentence about locus generators in front of the one
  sentence naming the setting and the limit.

## Where calling was put, and how the merge kept its oracles

**The builder calls each locus, and the merge module still knows nothing about calling.**
`build_region`'s locus walk moved into `build_region_handing_over`, which hands each surviving
locus to a sink and each refused locus's span to a vector; `build_region` is that function with
`Vec::push` for a sink. `merge_cohort_through_cache` split the same way into
`merge_cohort_handing_each_locus_over`. The run supplies the calling sink from `callers.rs`, so
`cohort_merge` imports nothing from `ng::calling`.

Three things that shape buys:

- **The ownership rule of spec §6.1 is written once.** Widening either end of that comparison
  loses a locus from the run with nothing to say so — about one in twenty at twenty-base
  building regions — so a second copy of the loop was never an option.
- **Every existing oracle of the merge is untouched.** `merge_cohort_through_cache` keeps its
  signature and its behaviour, so C2's differential and the merge's own 369 tests check the
  streaming driver and the collecting one at once.
- **`merge_cohort` stays**, as the merge's oracle rather than the run's path, and keeps none of
  the walk's tallies on purpose.

**The proof that the placement is free** is
`calling_inside_the_builder_gives_what_calling_after_the_merge_gives`: the same cohort called
both ways, refusing any difference. What it does **not** catch is a scratch that leaks between
loci — both sides walk the same loci in the same order, each reusing its own scratch in the same
pattern, so a missed `clear()` produces the same wrong answer on both sides. Spec §8's
calling-scratch trap is caught where the *order* differs, which is E2's.

## The time split D3 needs

Milestone E's shape is decided from where a real run's time goes, and until now the builder's
time was one undifferentiated number. `AFTER_ASSEMBLY_NANOS` times the sink — in a calling run,
the genotyping — so `Report::assembling_loci_ms()` is the assembling alone, and both now print
in the report's own breakdown. Read once per built locus, which costs two `clock_gettime` calls
a locus under `--features merge-timing` and nothing at all without it.

**One existing measurement changed.** The serial cached driver no longer adds to
`ORGANISE_NANOS`: it used to charge its per-region `Vec::extend` there and it no longer has one.
A serial breakdown now reads zero for "releasing loci in region order", which is the parallel
driver's own phase. The counter's own documentation says so.

## Verification

| check | result |
|---|---|
| `cargo test --lib` | 5,809 passed, 13 ignored (5,801 before this step) |
| `cargo test --lib ng::run` | 369 passed (361 before) |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo doc --no-deps` | 26 unresolved links, 23 redundant link targets — the standing baseline, unchanged |

## What the reviews changed

Three reviews ran in parallel over the step's diff, each in its own worktree: one mutation
testing, one on design fidelity against the spec and architecture, one on what a person sees.

### The design review found a regression this step introduced

**The merge's own oracle was walking with different generator settings from the run it checks.**
The new `locus_generator_settings` field was given the deliberately-unusual value at every
fixture, including the four that feed C2's differential — and that differential's oracle builds
its own generators with the shipped constants. The module header says the two share the
generator set and cannot not share it. Those four fixtures now use the shipped settings; the
unusual ones stayed where they belong, in the construction checks. It passed only because the
fixtures are three reads deep, where no depth cap of any plausible size binds.

### Five ways to be wrong with the suite green

The correctness review wrote 18 mutations and 6 lived. **Five are now dead, each re-injected and
re-checked after its test was written.** Four of the five produce wrong genotypes rather than a
crash:

- **The run's candidate selection, calling-loop settings and merge parameters could each be
  replaced by the shipped defaults** inside `call_cohort`, so an operator's own thresholds would
  be silently ignored. All three survived for one reason: the calling fixtures opened with the
  shipped values, so the run's settings and the defaults coincided and no substitution could
  change an answer. Three tests now open a run with a non-default setting and assert the answer
  moves — one called locus against two at a different minimum-alt floor, two alleles against one
  at a different candidate floor, one pass against several at a different pass cap.
- **The locus generator's settings could be thrown away** and every generator built from the
  constants, which is the whole of what the new field is for. `locus_generator_settings()` only
  proved that `open` had stored the value. The fixture that separates them puts the per-position
  cap at 1 against the shipped 8,000 on a three-read column, where `column_depth_truncations` is
  0 at one setting and non-zero at the other. **`positions_short_of_cap` cannot see it**, which
  is why it is not what is asserted: that counter answers only whether the read-hold ceiling
  cost coverage.
- **Every sample's walk tallies could be permuted across samples**, reversed on the way out of
  the cache or reversed at the join with the names. The test pinned the names in run order and
  every count it asserted was a property all three samples shared. The cohort now walks 3, 5 and
  7 reads, so the counts identify the sample they came from.

### One survivor could not be pinned, and the fixture reference is why

**A run's list of loci the merge declined to assemble can be replaced by an empty vector with
every test green.** Pinning it needs a locus wider than one reference base, which on real reads
means a deletion — and this module's fixture reference is a hundred `A`s on `chr1`, so every
deletion in it sits inside one homopolymer. Measured: three reads carrying a five-base deletion
at `chr1:20` produce **no cohort locus at all**, at the shipped width bound and at a bound of
three alike, so there is nothing for the width test to refuse. What would pin it is a fixture
reference with varied bases, which `build_fasta` does not build — it takes contig names and
lengths and fills with `A` — and which four modules share. Recorded at the point in the code
where it matters rather than changed under this step.

**And one survivor is not a defect.** The sink's own timing can be made to measure nothing and
no test notices, because nothing reads `Report::after_assembly_ms` and the counters are behind
`--features merge-timing` in any case. The honest statement is that the split is documented and
not verified; D3 is the first thing that will read it.

### Nine wrong claims in the new prose

The worst was mine about my own instrument: `positions_short_of_cap` at zero was written as
"the run's depth settings shaped no evidence", and it means only that the **read-hold ceiling**
cost no coverage. The two per-position caps are a different counter, and a run can have zero on
the first and millions on the second — the ceiling kept every read and the caps then declined to
score on all of them. A reader following that sentence would have concluded the opposite of the
truth.

The others: `RunError::LocusGeneratorSettings` still said it could never fire; `views` was
called "the one allocation per locus" when the call itself allocates a `LocusInference`'s
per-sample calls; an `#[expect]` reason counted "four things" in front of a list of five; a
clamp was justified by arithmetic that cannot occur (the sink runs strictly inside the builder's
stopwatch); the read-filter gap was blamed on a missing accessor rather than on the contig
boundary that drops the counts; `AlignmentInputs`'s own grouping rationale still described five
fields; "one called locus per surviving position" contradicted the rule that a locus is not a
position; and `CohortWalkTallies` opened by saying the walkers do not survive the merge, two
lines above the code that brings them back.

## What this step does not do

- **It does not yield one record at a time.** `call_cohort` returns every called locus at once,
  where spec §5.1 bounds a run at `callers in flight × one cohort locus` plus the frontier. What
  it *does* bound is the observations — each is dropped as soon as its genotypes exist — and the
  refused-span list accumulates for the whole run as well. The pool milestone is where the calls
  start being released singly, and it inherits a driver that no longer has to buffer the
  observations to get there.
- **It produces called loci, not VCF records.** The writer is coded and the format settled; D3
  is where a run's output becomes a record.
- **No repeat tract goes through it.** Both tract generator slots are unfilled, which is the
  plan's own scope decision, so a run over ground with tracts in it is short rather than wrong
  and the walk tallies say by how much.

## ⚑ One thing for the owner, before D3

**A locus where the allele cap has ruled every covering sample uncallable aborts the run with a
panic**, and `call_cohort` is the first thing that can reach it from real data
(`summarise_condition.rs`, *"this locus has nobody to call"*). Spec §4.1's ruling — selection
cuts an allele rather than refusing a locus, because most samples stay callable — explicitly
does not cover the case where none does.

It is unmeasured. The closest measurement argues it is rare: on HG002 at 300× no locus of 7,478
carried more than three alternatives against a cap of six, so the cap never bound there. At
tomato depth the share is inert and the bar is a count of two, so alternatives accumulate across
63 accessions and the cap can bind — though 63 samples also means many non-covering samples,
which stay callable. D3 is the first run that can answer it with a number.
