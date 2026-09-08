# Research plan — what the cohort merge's parallel driver actually costs, and what would fix it

**Date:** 2026-08-28. **Status:** plan, no work done. **Owner's ask:** the parallel merge's
speed-up is poor; find out why and try the candidate fixes.

**This is a measurement job first and a coding job second.** Nothing here is to be adopted on the
strength of an argument. The deliverable is a written finding — including the changes that were
built and refuted, which are as valuable as the ones that worked.

---

## 1. The question

`merge_cohort_in_parallel`
([`src/ng/run/cohort_merge/parallel.rs`](../../../../src/ng/run/cohort_merge/parallel.rs)) does the
same work as the single-threaded merge and produces byte-identical output, but gives back much
less speed than the thread count suggests. **Why, and what would recover it?**

**The size is known and so is the leading cause.** Eight threads give **1.4×** on 63 tomato
accessions, and a narrow enough region makes threading *slower than not threading*
(`src/ng/run/cohort_merge/mod.rs:545-571`, and §2 below). The same comment says why: the merge works
in **rounds** — the organiser draws every sample's reader forward while no builder runs, then that
round's builders run together, then nothing is released until the slowest of them finishes — and
*"the cost a narrow region adds falls on the organiser, which no thread but one ever runs."*

**So the question is not "what is the suspect" but "how much of the 1.4× gap is the organiser, and
what else is in there".** Two other costs are on the record and neither is the barrier: a fixed
per-region cost that walks the whole cohort, and a 39% allocator share (§2). Step 1 sizes all three
before anything is changed.

---

## 2. What is already known — read this before measuring anything

**Two probes already exist. Do not write a third until you have run both.**

- [`examples/ng_cohort_merge_parallel_cost.rs`](../../../../examples/ng_cohort_merge_parallel_cost.rs)
  — times the oracle, the cached serial driver and the parallel driver over fabricated ground, and
  sweeps the building region's width and the cohort size.
- [`examples/ng_cohort_merge_real_cost.rs`](../../../../examples/ng_cohort_merge_real_cost.rs) —
  runs the merge on observations real reads produced, and times the locus generator beside it so
  the merge's share of a run can be seen.

**Their headers carry findings that are written nowhere else**, which is one reason this plan
exists. What they say:

- **The building region's width has already been swept on real reads, and the answer is recorded in
  the code rather than in any document** (`src/ng/run/cohort_merge/mod.rs:545-571`). On 63 tomato
  accessions over 100 kb of SL4.0: **one thread barely notices the width** — 656 ms at 20 bases,
  615 at 100, 616 at 200, 624 at 500, 636 at 1,000 — while **eight threads notice a great deal**.
  At 16 samples, eight threads take **173 ms at 20 bases against 130 ms on one thread**, so
  threading a 20-base merge is *slower than not threading it*, and **93 ms at 200 bases, which is
  1.4× one thread**. The eight-thread optimum is 100–200 bases at both cohort sizes, and 200 is
  the shipped default.

- **That measurement already names the cause**, in the same comment: *"the cost a narrow region
  adds falls on the organiser, which no thread but one ever runs."* So the serial phase between
  rounds is not a hypothesis about where the time goes — it is the recorded explanation of a
  measured effect, and 1.4× on eight threads is the size of what is left on the table.
- **The merge pays a fixed cost per building region per sample** — a cover, an eviction, a window,
  and a builder's setup, each of which walks the whole cohort. That cost is serial-ish work that
  grows with cohort size and has nothing to do with the barrier.
- **The merge's Linux profile is 39% allocator.** If that still holds, allocation may be a larger
  lever than scheduling.
- **Real observations are about one record per base per sample** on the tomato benchmark's 63
  accessions — a hundred times denser than the "sparse" fabricated fixture, and denser than the
  "dense" one.
- **Nobody has ever timed the merge against the rest of a run.** It may not be worth optimising at
  all; the real-cost probe prints the generator's time beside it for exactly this reason.

**Measurement hygiene, learned the hard way and already in those headers:** report the median of
five repeats with the spread; this machine drifted 14% across two runs of one unchanged binary an
hour apart, and swung 30% between runs at 3,000 samples. Compare two builds by running them
alternately in one sitting, never one after the other. Build the fixture outside the timed span —
charging one driver for a clone the other never did made the cached path look two to three times
worse than it is.

---

## 3. Step 1 — where the time goes (measure, change nothing)

**Produce a breakdown of the parallel merge's wall time into five parts**, on real observations,
at a few cohort sizes:

| part | what it is |
|---|---|
| the organiser's serial phase | drawing the readers forward between rounds, while no builder runs |
| the barrier | builders finished and waiting for the slowest in their round |
| builder work | the merge proper |
| the per-region fixed cost | cover, evict, window, builder setup — each a walk over the cohort |
| the allocator | what fraction of all of it is allocation and free |

**This is the step that decides everything after it.** The organiser's serial phase is the leading
suspect on the record, but its *share* has never been measured against the other two, and the
ordering of §4 and §5 depends entirely on that split. If the allocator turns out to dominate, §4b is
the whole job and the structural work in §5 is wasted.

**How.** Instrument the parallel driver behind a feature or a compile-time switch — timing that
ships in the hot path is itself a change to what is being measured. The allocator share wants a
profiler rather than counters (`samply` on macOS against the host binary; on the Linux box, see
the profiling notes in `CLAUDE.md` — `perf_event_paranoid` has to be lowered on the *host*, and no
container flag substitutes).

**Data.** `benchmarks/tomato1/crams/` holds 63 sliced CRAMs of about 11 MB each;
`benchmarks/tomato1/regions.bed` and `regions_n160_200kb.bed` are the interval sets. The reference
lives under `$HOME/genomes` and is mounted read-only in the dev container.
`ng_cohort_merge_real_cost.rs` already takes `<reference.fa> <cram-dir> <regions.bed>` and honours
`NG_REAL_SAMPLES` and `NG_REAL_REGIONS` — use those to sweep the cohort size rather than making
new fixtures.

**Also report the merge's share of the two stages together** — generator plus merge — because if
the merge is a tenth of the run, the ceiling on all of this is 10%.

---

## 4. Step 2 — the two cheap moves

Both are candidates the profile of §3 may rule out. Do not do them because they are listed.

### 4a. Extend the width sweep — do not repeat it

**The sweep exists on real reads and §2 gives its numbers; repeating it is waste.** What it does
not cover, and what is worth an afternoon:

- **more than eight threads.** Everything recorded is one thread against eight. If the organiser's
  serial phase is the cap, the curve should flatten and then turn down, and where it turns is the
  number that says how much §5's structural work can possibly buy.
- **the far end of the cohort range.** The recorded run is 16 and 63 samples. The organiser's
  per-region work walks the whole cohort, so its share should *grow* with sample count — meaning
  the 1.4× at 63 samples may be the optimistic figure, not the pessimistic one. Sweep toward a
  thousand with `NG_REAL_SAMPLES`.

**Report peak resident beside wall time**: the cache holds `regions in flight × width` bases plus
the tail reaching past it, which is the whole reason the regions are short.

### 4b. Take the allocations out of the per-region path

If the 39% allocator share holds, this is the largest lever and **it does not touch the schedule at
all** — it makes the serial phase, the builders and the barrier all cheaper at once, and it is
orthogonal to §5, so it can be adopted whatever the structural work concludes.

The shape: buffers that live for the run and are cleared and refilled per region, rather than
allocated per region per sample — the load / use / clear / reload pattern this project already
prefers. The per-region fixed cost the probe names (cover, evict, window, builder setup, each a
walk over the cohort) is where to look first.

---

## 5. Step 3 — the structural changes, if the profile justifies them

### 5a. Overlap the reader advance with the building

The organiser draws the readers forward while no builder runs, and that phase is input and
decompression — the work most worth hiding behind computation. Fill the **next** round's ground
into a second buffer while the current round's builders work, then swap.

**Cost:** two rounds of observations resident instead of one, so roughly double the cache's peak.
Measure it, do not assume it.

**Correctness needs no new argument:** release is still by region index, so the output order cannot
change. §7's oracle still has to pass, and if it does not, the bug is in the buffering.

### 5b. Drop the rounds for a sliding window

**The barrier is not there for ordering.** The organiser already releases along a gapless run of
region indexes, which supports builders finishing in any order at all
([`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs)). It is there because builders
read the observation cache while the organiser writes it, so the two cannot overlap.

Make the cache a window: extended ahead of the claim frontier, evicted behind the release frontier,
with a published "covered to here" mark that builders never read past. A builder that finishes
takes the next unclaimed region immediately, whatever its neighbours are doing. No rounds, no
barrier.

**This is the largest of the changes and it subsumes 5a**, so it is the one to reach for only when
the earlier steps have said the barrier is what costs. It also changes the memory story — the
window's bound is no longer `regions in flight × width` and has to be re-derived and re-measured.

---

## 6. How to work

**Take your own worktree and branch.** From the repository root:

```
git worktree add ../pop_var_caller-ng-merge-parallel-cost -b ng-merge-parallel-cost
```

Work only inside that worktree; do not touch the main checkout. The convention is
`../pop_var_caller-<feature>` with a plain `git worktree add`.

**Building and running.** On macOS, everything that compiles goes through `./scripts/dev.sh` by
absolute path. On the Linux box `rick` there is no container runtime and `cargo` is run directly —
and the two put their binaries in different places (`target-container/release/` against
`target/release/`), so a script that looks for a built binary must check both and take the newer.
`CLAUDE.md` has the details, including why profilers need a sysctl lowered on the host.

**Scratch files go in the repository's own `tmp/`**, never in `/tmp` and never in a harness
scratch directory: a path outside the project mount is invisible inside the dev container.

**Do not edit `src/ssr/`, `src/pileup/`, `src/psp/`, `src/var_calling/` or `src/vcf/`.** Production
is frozen. This work is confined to `src/ng/run/cohort_merge/` and `examples/`.

---

## 7. What must stay true

**The one hard constraint: the output cannot change.** The parallel merge's whole claim is that it
gives the same answer as the single-threaded one, byte for byte, at any building-region width and
any number of regions in flight. Every change here must keep the existing suite green, and any
change to the schedule must be checked against the serial oracle explicitly, not only against the
suite.

**A change measured only on the fabricated fixture is not adopted.** That fixture is a hundred
times sparser than real ground and puts every sample's records at the same positions, which real
samples do not do. Fabricated numbers are for bracketing the shape of a cost across a wide density
range; real numbers decide.

**Report what was refuted.** A change that was built, measured and set aside is a result. The
archived experiments folder exists for exactly this
([`experiments/README.md`](experiments/README.md)) and says why: the numbers that were rejected and
the code that was deliberately not adopted are what a conclusion cannot carry on its own.

---

## 8. What "done" looks like

1. **A finding document in this folder**, `cohort_merge_parallel_cost_<date>.md`, answering: what
   the merge's time is made of, what share of a run it is, which of the four candidates were tried,
   what each was worth, and which are recommended. Numbers with their cohort size, their density,
   their machine and their allocator — a figure measured at 63 samples on sliced CRAMs is a fact
   about that corner, not a property of the merge.
2. **The adopted changes as commits on the branch**, each keeping the byte-identity oracle green,
   with the throwaway measurement code kept separate from anything proposed for merge.
3. **`doc/devel/ng/spec/run_streaming.md` §11, question 7 answered** — that question names these
   three ideas and is where the conclusion is owed. If the answer is "the merge is a small share of
   the run and none of this is worth doing", say that; it is the most useful outcome on offer and
   it closes the question for good.

---

## 9. What is deliberately not in scope

- **Parallelising anything other than the merge.** Whether the decode or the calling loop deserves
  a pool is the same open question's other half and needs the psp format, which does not exist yet.
- **Changing what the merge computes.** The keep rule, the locus span bound, the allele table — all
  fixed. This is about how the same answer is produced.
- **Adopting a default.** Recommending a building-region width is in scope; setting it as the
  shipped default is the owner's, and `run_streaming.md` currently leaves every concurrency default
  unset on purpose until the psp format's costs are known too.
