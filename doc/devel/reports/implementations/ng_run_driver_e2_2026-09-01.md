# E2 — concurrency invariance: the same VCF at every thread count

**Date:** 2026-09-01. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone E step E2, and **Checkpoint E**. **Design:** [`run_streaming.md`](../../ng/spec/run_streaming.md)
§12 oracle 2, §8 (the calling-scratch trap). **Module:** `src/ng/run/callers.rs`, test module
`records_handed_over_as_the_run_finishes_them`.

---

## What the oracle pins, and what it cannot pin yet

**Spec §12.2 asks for the same VCF at any number of callers in flight; under E1's shape that
means any rayon thread count**, the cover's sweep being the only thing a calling run threads.
The plan wrote E2 around the calling-scratch trap — a missed `clear()` that shows only when a
pool reorders the loci a scratch sees — and **that trap cannot fire yet**: no built arrangement
reorders the loci, because assembly and genotyping stay on the merge thread in genome order.
The oracle is built anyway, so whatever first threads the calling is caught by a test that
predates it rather than owed one.

## What landed — three tests

**`the_record_path_is_byte_identical_at_every_thread_count`.** A five-sample cohort whose loci
differ in kind — a SNP two samples share (three and four reads, so the samples are numerically
distinguishable), a sample-private SNP, a two-base insertion, a locus that is called and
establishes no variant, and a sample with only reference reads — over ground with a twelve-base
`AT` repeat tract interleaved between two ordinary stretches, which is the plan's "ordinary
sites and repeat tracts interleaved" as far as the built caller can differ (the tract routes to
the unfilled slot and is charged to *not built yet*; the test pins that the routing and the
accounting are thread-count-invariant too). The cohort's VCF is written through the real
`VcfWriter` at pools of 1 (the serial-sweep fallback branch), 2, 4, 8 and 16, three repetitions
per parallel pool, and the files are compared **byte for byte** — along with every count a run
report is built from: records written, called-but-not-written, the two refusal lists — **both
empty on this fixture, so those two comparisons carry no weight and are kept as a guard against
a future fixture that fills them** — and each
sample's walk tallies.

**`the_mixed_cohorts_records_describe_the_serial_callers_loci`** ties the sweep back to the
oracle the plan names: `call_cohort`, which never touches the parallel cover, calls exactly the
loci the record path writes plus the called-but-not-written count. It runs the record path in a
named pool of eight and asserts it, so the comparison cannot degrade to the serial cover
against itself on a one-CPU runner.

**`a_cohort_of_one_sample_is_byte_identical_at_every_thread_count`** takes the same sweep down
to the hardest end of the range this caller commits to. One sample means the sweep's reduction
folds a single value, so anything that only works because a second sample happened to widen the
reach shows here and nowhere else in this module.

## What the width sweep is, and what this test cannot catch

The byte test runs everything at **two building-region widths** — the shipped default and seven
bases — and asserts the bytes are identical across the widths as well as across the pools. At
seven bases the fixture's 40-base first stretch is divided instead of taken whole, so the run
covers, evicts and builds its ground in several steps while answering the same file. That is
what the second width buys: **not a stronger oracle, a different division of the ground**.

**⚑ Two corrections about that width, both from this step's review.** The shipped default is
**500 reference bases** (`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`,
`src/ng/run/cohort_merge/mod.rs:592`), not the 200 an earlier draft of this report and
`cohort_merge.md` §6.4 both say — the 200→500 change is the research note's adopted candidate
and the spec is the owner's to bring into line. And **nothing asserts that the width a caller
is given reaches the merge**: replacing it with the constant default at either calling entry
point leaves all 437 tests green. Width is output-invariant by design, so no comparison of
outputs can close that; it would take counting the covers under `--features merge-timing`.
Recorded rather than claimed.

**⚑ It does not make the test able to see the cover's fixpoint, and an earlier draft of this
report said it did.** Three mutations to `cover_in_parallel`
(`src/ng/run/cohort_merge/observation_cache.rs`), each applied alone and reverted before the
next, run in the container against this step's two tests:

| mutation | effect on the sweep | these two tests | `ng::run` as a whole |
|---|---|---|---|
| drop the last sample from `par_iter_mut` | a sample is never drawn forward | **both fail** | fails |
| `max` → `min` in the `try_reduce` | a cover stops at the least reach any sample grew to | 2 passed | 431 passed, **6 failed** |
| `break` after the first iteration | the fixpoint never iterates | 2 passed | (not run whole) |

```
$ ./scripts/dev.sh bash -c 'cargo test --lib -- \
    the_record_path_is_byte_identical_at_every_thread_count \
    the_mixed_cohorts_records_describe_the_serial_callers_loci'
test result: ok. 2 passed; 0 failed; 5919 filtered out          # with min-for-max
test result: FAILED. 0 passed; 2 failed; 5919 filtered out      # with a sample dropped
```

**So the oracle has real power over *who* the sweep draws and none over *how far* it keeps
drawing.** The six tests that do catch the fixpoint mutations are
`the_parallel_cover_gives_the_serial_drivers_answer`,
`the_cache_changes_nothing_at_every_building_region_width`,
`the_two_drivers_agree_on_random_layouts` and the three parallel-driver comparisons:

```
$ ./scripts/dev.sh bash -c 'cargo test --lib ng::run'
test result: FAILED. 431 passed; 6 failed                       # with min-for-max
```

**The reason is the fixture reference, and it is not fixable here.** These tests call from real
BAM files against the shared fixture reference, which is a hundred identical bases. A cover has
to follow a chain only when one sample's observation reaches past a building region into where
another sample's begins, and on that reference nothing reaches: two substitutions never share a
base, an insertion's reference span is its anchor alone, and a deletion left-aligns off the
record (measured at D1). Two samples departing at adjacent positions were tried during this
step and closed as **two separate loci**, not one. With no chain, there is nothing for a second
sweep to find, so a cover that stops early loses nothing on this ground — which is why the two
fixpoint mutations are invisible here and a dropped sample is not. `cohort_merge`'s own
fixtures are minted in memory and can hold a 26-base observation, which is why the chain claim
belongs there.

So the division of labour is: **the cover's schedule-independence is pinned a layer down, on a
fixture that can express a chain; what this step adds is the end-to-end tie** — the whole run,
through the real writer, to the bytes a person diffs. `the_mixed_cohorts_records_describe_the_serial_callers_loci`
now asserts the one-position width, so the limitation is checked rather than assumed: if a
later change ever does produce a wider locus here, the test fails and this section has to be
rewritten instead of quietly becoming false.

## How this step was reviewed, and what the reviews changed

The fixture was written against the E1 determinism review's probe suite, which ran the same
end-to-end shape at pools of 1, 2, 4 and 8 and found the bytes identical. That inheritance is
not a review of this step, so the step went through the per-step fan-out — two agents in
isolated worktrees, each re-pointed at ee7124f0 with the diff handed over as a patch.

**The determinism attack could not make the output differ, and what it threw is the evidence
behind that.** Roughly 19,000 comparisons: 15,288 parallel cover traces against 392 serial ones
over 1–40 samples × 4 chain shapes × 14 building-region widths × 13 pool sizes; 1,680 merge
runs against the serial oracle over 120 random layouts; 810 whole runs to VCF bytes over five
cohorts (including one sample, and ten samples of pure allele ties where an `ahash` iteration
order would show) × 6 widths × 13 pools; and repeats of the cover and end-to-end sweeps under
16 busy threads. It compared not only the bytes but `covered_to`, the held count, each sample's
recycle-list length and the full `Debug` of every held observation after every building region.
All identical.

**The correctness pass ran 18 mutations and found one real hole on the path the command
runs**, which is now closed. Reversing the sample-name list before it is paired with the walkers
in `call_cohort_handing_each_record_over` passed all 437 tests; the same reversal in
`call_cohort` was killed by two. The guard existed only on the oracle path. A run report would
have carried one sample's read-drop rates under another sample's name — a wrong report rather
than a crash, and `CohortWalkTallies::of`'s own length check cannot see it because both lists
are the right length. `the_record_path_is_byte_identical_at_every_thread_count` now asserts
each sample's admitted reads under its own name against absolute numbers (3, 4, 3, 3 and 6, so
no two samples are interchangeable), and re-running the mutation fails it.

Four further changes came out of the two reviews:

- **The per-sample comparison was comparing three numbers that are identical for all five
  samples**, so only the name column varied and it never varies with thread count. It now
  compares `LocusCounts`, the read-filter tallies and `PileupGeneratorCounts` **whole**, so a
  count added to any of them later joins the comparison instead of dropping out of it.
- **The oracle-tie test was not pinned to more than one thread.** It ran under the ambient
  pool, so on a one-CPU runner or under `RAYON_NUM_THREADS=1` it would have compared the serial
  cover with itself and still passed. It now runs in a named pool of eight and asserts it.
- **A cohort of one sample is now swept across thread counts** — the shape the project's own
  range commitment makes the hardest and the one where a sweep over samples has least to do.
- **Two doc comments contradicted each other** about the width sweep, one of them carrying the
  retracted mutation claim. Both are rewritten.

Three findings were recorded rather than fixed, because no comparison of outputs can close
them: the building-region width need not reach either calling entry point for the tests to pass
(width is output-invariant by design); both refusal lists are empty on this fixture; and at the
shipped width the parallel sweep runs once per run, since every read of the fixture lies in the
first of the three typed regions.

## What was measured

In the development container at this step's working tree, all commands quoted in full.

| check | result |
|---|---|
| `cargo test --lib` | **5,908 passed, 0 failed, 14 ignored** — 5,905 at ee7124f0 + the 3 added here |
| `cargo test --lib ng::run` | **438 passed** — 435 at ee7124f0 + the 3 added here |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 (checked by exit code, not by reading a tailed log) |
| `cargo doc --no-deps` | 26 `error: unresolved link`, 23 `warning: redundant explicit link target` — the standing baseline, added to neither |

The baseline row of that table was re-measured on a clean ee7124f0 at the start of this
session rather than quoted, and matched: 5,905 passed / 14 ignored, `fmt` exit 0.

## Two things this session found that were not this step's work

**The development container's incremental-compilation cache was corrupt, and it forged a
determinism failure.** Before any of this step's code was touched, `cargo test --lib ng::run`
at ee7124f0 returned **429 passed, 6 failed** — the six merge tests above, with the parallel
side assembling fewer of a sample's records into a locus. The same cache then crashed rustc
with an internal error in `join_codegen`, produced a link failure against undefined
`anon.*.llvm.*` symbols, and produced a binary the container could not exec. Deleting
`target-container/debug/incremental` (2.6 GB) ended all four symptoms. The six tests have since
passed in **457 consecutive runs**, 331 of them with all eight container CPUs saturated to force
contention. The signature of that false failure is *identical* to the signature of the `min`
mutation above, which is what identifies it as a miscompiled binary rather than a race.

**The earlier draft of this report attributed its own lost work to the owner's commit.** It did
not: the module was reverted by this session's own `git checkout` while establishing a clean
baseline, the uncommitted work having been mistaken for a dirty tree. Recorded because the
report stated a cause it had not checked.
