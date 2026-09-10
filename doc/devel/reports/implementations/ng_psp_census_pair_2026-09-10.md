# The census lives inside the psp — the whole plan

**Date:** 2026-09-10
**Plan:** [psp_census_pair.md](../../ng/impl_plan/psp_census_pair.md), Milestones A to E
**Spec:** [psp_census_pair.md](../../ng/spec/psp_census_pair.md)
**Branch:** `census-vs-psp-perf`, from `9fdc427f` to Checkpoint E — plus the two measuring
programs under `examples/`, which landed on the branch before that commit and which Milestone E
checks rather than lands

This is the plan's own report, written at its end. Each milestone has its own beside it —
[A](ng_psp_census_pair_milestone_a_2026-09-09.md), [B](ng_psp_census_pair_milestone_b_2026-09-09.md),
[C](ng_psp_census_pair_milestone_c_2026-09-10.md), [D](ng_psp_census_pair_milestone_d_2026-09-10.md),
[E](ng_psp_census_pair_milestone_e_2026-09-10.md) — and a review report each. What is here is what
a reader needs who was not following along: what changed for a person running the caller, what it
cost, and what was decided along the way that the spec did not already settle.

---

## What a person sees now

**A sample is one file.** Before this plan a walk left two: the psp, holding what that sample's
reads showed at every position of the ground, and beside it a `<sample>.census` — the much smaller
object a parameters fit reads. Now the census is sealed into the psp's tail, and there is no
second file.

The pipeline is three commands and a repair:

    generate-psps        alignments        ->  <sample>.psp, census sealed inside
    estimate-parameters  psps              ->  cohort.parameters.toml, or a refusal naming every stale sample
    call-from-psps       psps + that file  ->  the VCF
    regenerate-census    psps              ->  each psp's trailer replaced, for the samples the fit refused

**Four commands became three and a repair.** The walk already built the census in the pass it
was making anyway; what `generate-census` did was build a *second* copy, from the stored psp, into
the file beside it. Every run made both, and an end-to-end run compared them — that agreement
between two producers is the oracle this plan kept and strengthened. What the plan removed is the
second copy, not the work: there is one census a sample now, in the psp, and the command that
built the other is `regenerate-census`, which is not a step of a run at all. It is what you run
when the fit refuses a psp, and it rewrites only that file's tail.

**No command takes a repeat criterion any more.** What counts as a repeat tract, what ground was
walked, which reference and catalog were used — a psp records all of it, and the fit and the repair
read it from there. Before, a person retyped those settings at step 2, and a mistyped one produced
a plausible census fitted under settings they thought they had not changed.

## The three things it makes impossible

**A census cannot come apart from its psp.** It used to be a separate file that could be copied
alone, deleted, or paired with a psp it was not built from. The census carried a digest of its
psp's header and that psp's record count so the mismatch could be caught; both are now deleted,
because a census that *is* the psp's tail has nothing to be paired wrongly with.

**A cohort cannot be fitted or called across psps that disagree on what they cover.** Every
command that opens a cohort compares each psp's header against the first's — the analysed regions,
the repeat catalog, the repeat-tract criteria — and refuses the set, naming the sample and the
field. `call-from-psps` gets that for free. **The read filters each walk applied are not among
them**: they are recorded in every psp's header and the calling run's report names the ones that
disagree, but nothing refuses over them.

**A stale census cannot be read as though it were fresh.** Twelve settings travel with each
census — nine of them saying which positions were chosen, three saying in what units the evidence
at those positions was written down — and the fit compares all twelve against what its own run
records under before it fits anything; a census that does not match is refused by name, with the
command to run. Before, the fit compared the digest of the kept positions alone, which let a
census built under other settings through whenever both happened to keep the same positions.

**Two of those refusals come at different moments, and the cheap one comes first.** Whether a psp
carries a census at all, and whether it is of a format this build reads, is ten bytes and a seek —
so the fit judges every psp that way **before it opens the reference**, and a cohort of sixty stale
psps is refused in the time it takes to open sixty files. Whether a census was recorded under other
settings cannot be known that cheaply: the settings are a digest over a selection of positions, and
the selection has to be rebuilt from the reference to have anything to compare. That refusal
therefore costs the reference read — a few seconds to twenty, on the references this caller is
used on.

## What it cost, measured

**Nothing a run produces moved.** On the first six tomato accessions of
`benchmarks/tomato1/crams` by name, over the first two 100 kb intervals of its `regions.bed`, at
about three reads a position, the whole pipeline gives the same census total to the byte —
**1,545,479** — the same parameters file size — **38,124 bytes** — and the same VCFs: **2,275
records with the compiled-in defaults against 2,082 with the fitted numbers**, 196 called only by
the defaults and 3 only by the fit, 87 of the 2,079 both called differing in at least one genotype,
**113 genotypes of 12,474 compared**.

**Same as what**: the run of 2026-09-09 recorded in
[`ng_psp_census_pair_milestone_a_2026-09-09.md`](ng_psp_census_pair_milestone_a_2026-09-09.md),
which is the comparison that matters, because that run's fit read census *files* and this one's
reads psp trailers.

**A psp whose census is rebuilt is the walked file byte for byte, whole** — header, blocks, index,
trailer, footer — on all six accessions on real reads, and on the fixture cohort in a test.

**Building the census during the walk is the cheaper of the two routes** — 1.28 s against 1.40 s
over six accessions on 200 kb, measured at the fit stage's Milestone B
([`ng_fit_stage_b_2026-09-05.md`](ng_fit_stage_b_2026-09-05.md)) and quoted here rather than
re-run. That was already the reason the walk builds one; this plan did not change which route a
run takes.

**One figure did move and it is not the caller.** The psp total came back 36 bytes larger than the
recorded run's, over six samples. A psp's header stores the run's command line verbatim, so the
total depends on how long the arguments were — measured: walking one accession into a path five
characters longer gave a psp five bytes larger, 1,787,534 against 1,787,539. Thirty-six bytes over
six psps that share one command line is six characters, and the earlier run's report does not
record the path it used, so which token was shorter cannot be said. The census inside those psps is
identical to the byte.

## What was decided during the build that the spec had not settled

**A census recorded under other settings is refused rather than rebuilt** — the owner's ruling
of 2026-09-10, which added a step to the plan. The alternative was to rebuild silently, and the
argument against it is that a rebuild costs a full pass over every psp: a person who did not ask
for one should be told, not billed.

**The refusal names the setting that differs, and says which of two fixes it calls for.** A wrong
reference or catalog is fixed by rerunning with the files the psps were walked against, which is
free; anything else needs `regenerate-census`, which costs a pass. A message that could not tell
them apart would leave the reader to guess.

**A census's format did not change when the identity was deleted.** The flag byte that said whether
a census named its psp is still written, always zero. Removing it would have moved every later
field, cost a format version, and made every psp already on disk unreadable — to save one byte a
sample. A census that *does* set the byte is refused: this build no longer has anything to compare
its naming against.

## What is left open

Listed at Checkpoint E in the handoff, and none of it blocks a run:

- **No golden census on disk.** Every round-trip test writes with this build and reads with this
  build, so a coordinated change to both sides of the format would pass every test while making
  every psp ever written unreadable.
- **Three ways of getting a psp's header digest have no reader outside one test** —
  `psp::header_digest` and `psp::header_and_its_digest` in the psp module, and
  `WriteStats::header_digest` on what the writer reports — their only consumer having been the
  deleted identity. Removing them is tidying rather than a saving: the walk hashes a header once a
  sample, against a walk that runs for minutes.
- **Two documents outside this plan describe the sidecar in the present tense**:
  `doc/devel/ng/arch/run_streaming.md`, which cites two deleted items by line, and
  `doc/devel/ng/spec/run_streaming.md`, in four places — §2, §5.2, §6.1 and §13's reuse map.
- **A psp that will not open stops the cohort at the first file**, so spec §4.1's *every sample is
  examined* holds for stale censuses and not for unopenable ones.
- **`cargo doc` is red** with 40 unresolved intra-doc links, none of them this branch's, and is not
  in the gate.

## The gate, unchanged from where the branch started

Four gates, judged as a set rather than against green, because the tree was red in four ways before
this branch and still is — none of them this plan's:

| gate | at the branch point | at Checkpoint E |
|---|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | 6,738 lib tests pass; one integration target red | 6,737 pass, the same target red |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 11 errors of 5 kinds in 6 files | the same 11, in the same 6 files |
| `cargo check --all-targets --keep-going` | 4 examples do not compile | the same 4 |
| `cargo fmt --check` | 4 files | the same 4 |

**The two fours are different sets**, overlapping in one file: four examples that do not compile,
and four files — three examples and one source file — that are not formatted.

The one lib test fewer is arithmetic: three deleted with the type they exercised, two added by
Milestone E's review.
