# E1 — the calling run's cover goes parallel across samples

**Date:** 2026-09-01. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone E step E1 — **reshaped by the measurement the milestone's own deferral note demanded**
(owner's rulings of 2026-08-31 and 2026-09-01: what E becomes is decided from where the time
goes). **Design:** [`run_streaming.md`](../../ng/spec/run_streaming.md) §3.5 (provisional until
measured, by its own first sentence), §12.2; [`arch/run_streaming.md`](../../ng/arch/run_streaming.md)
§2, §3.1, §8 (amended by this step).

---

## The measurement, first — the whole cohort, both grounds

The plan's instruction was to measure at a cohort size that decides the milestone before
building it. Measured in the development container, release build, `--features merge-timing`,
`examples/ng_call_cohort_end_to_end.rs` (`NG_SAMPLES=63`, `NG_REGIONS=2` and `NG_REGIONS=80`),
2026-09-01, on the session's pre-E1 tree (the checkpoint measurement runs both arms on the
finished tree, so the serial arm is re-measured there):

| ground | loci called | `call_cohort` | drawing the readers | evicting | assembling | genotyping |
|---|---|---|---|---|---|---|
| 63 samples, 200 kb (D3's own ground) | 23,450 | 20.55 s | 18,105 ms (88.1%) | 172 ms (0.8%) | 1,135 ms (5.5%) | 1,091 ms (5.3%) |
| 63 samples, the whole 8 Mb benchmark | 1,069,772 | 723.2 s | 635,310 ms (87.9%) | 5,727 ms (0.8%) | 39,389 ms (5.4%) | 42,711 ms (5.9%) |

Both runs on one thread: 23.2 s of user CPU against 23.7 s elapsed on the first, 12m05s against
12m13s on the second (`time`, in the container). The split D3 measured at 3–24 samples holds at
63 and holds across a 40-fold change of ground.

**What that decides.** The plan's E1 — a pool of workers that genotype while the merge stays on
one thread — parallelises the genotyping column: 5.3–5.9% of a run. Its ceiling, with infinitely
many workers, is 1.06×. Switching the merge's region batching on instead parallelises assembling
and genotyping together: 10.8–11.3%, ceiling 1.13×. **Neither is worth its machinery at any
cohort size anyone can measure today.** The column worth parallelising is the first one, and
`ObservationCache::cover_in_parallel` — which sweeps the cohort's samples concurrently and
reaches the same fixpoint by a different schedule — already existed, reachable only by the
merge's parallel driver. So E's first question ("may a calling run use the parallel cover") is
answered **yes, and that is the whole of what E1 builds**; its second (which genotyping
arrangement) is answered **neither, at this cohort size** — revisit when a run big enough to
move genotyping's share is measured, which D3's report bounds at a tenth to a third of a
thousand-sample run across three defensible models.

## What landed

**The run's record path covers each building region with the samples swept concurrently;
everything else stays on the merge thread, in genome order, unchanged.**

1. **`merge_cohort_handing_each_locus_over_covering_samples_in_parallel`**
   (`src/ng/run/cohort_merge/serial.rs`) — the handing-over driver with
   `cover_in_parallel` in place of `cover`, taking the serial sweep on a one-thread pool the
   same way the parallel merge asks `rayon::current_num_threads()` before handing eviction to
   workers. The two public drivers share one private body
   (`merge_handing_each_locus_over_with`), so they cannot drift on anything but the sweep
   schedule.
2. **`AlignedFilesVariantCaller::call_cohort_handing_each_record_over` uses it** — the path
   `call-from-alignments` drives. **`call_cohort` keeps the serial cover on purpose**: it is
   the oracle the parallel-covered path is compared against, and one driver has to keep the
   schedule the fixtures were reasoned about under.
3. **The walk stack became `Send`**, which is what the sweep needed. Two of the three sites'
   own documentation had reserved room for exactly this change; the third had not, and the
   distinction is kept because it is doing work:
   - `GeneratorSlot::Generator` is `Box<dyn LocusGenerator<S> + Send>` — its doc named the
     exact condition ("if a `GeneratorSet` is ever moved onto a producer thread … this becomes
     `dyn LocusGenerator<S> + Send`");
   - the pileup generator's read-preparation cell is `Arc<Mutex<ReadPreparation<P>>>` where it
     was `Rc<RefCell<_>>` — its doc recorded the `Rc` as the remaining `Send` blocker;
   - `WindowedRefSeq`'s resident window is behind a `Mutex` where it was a `RefCell`, so the
     type is `Sync` and the `Arc<WindowedRefSeq>` handles the generator keeps are `Send`.
     **This one's doc had anticipated nothing** — it stated "`Send` but not `Sync` (per-worker
     ownership, like the production fetchers)" as the design. What makes the swap safe is the
     per-worker ownership that sentence records, not a reserved door: the lock is uncontended
     because nothing shares an accessor, before or after.
   `a_run_walker_can_cross_a_thread` is the compile-time proof, and it covers every filling of
   the generator slots at once because the box carries the bound.
4. **The probe grew the comparison arm**: `NG_COVER=serial|parallel` on
   `examples/ng_call_cohort_end_to_end.rs` selects the oracle or the run's own path, so the
   parallel cover's worth is measured by the same tool that produced the table above.

## Why the locks cannot change what a run answers, and what they cost

**Ownership is unchanged.** One walker is used by one thread at a time — the sweep moves it,
nothing shares it — so every lock introduced is uncontended. The trap the spec's §8 records
(`run_streaming.md`, "wrapping a shared one in a lock serialises the walk's hottest path") is about
*sharing* one accessor between workers, which nothing here does: each walker still owns its
three accessors.

**The schedule cannot change the answer**, and that argument is the cover's own, not this
step's: `cover_in_parallel`'s documentation carries the fixpoint argument (Jacobi against
Gauss-Seidel, same least fixpoint, same held window), and the parallel merge's entire oracle
battery has rested on it since Milestone E of the merge plan. What this step adds on top is
pinned by `the_parallel_cover_gives_the_serial_drivers_answer`: the handing-over driver against
the undivided oracle, five widths × pools of 1, 2 and 8, on the fixture carrying a
boundary-straddling locus and a width-refused span.

**One failure shape is less determined than before**: two samples failing in the same sweep race
for which error returns. The serial sweep always reports the first in sample order; the parallel
merge has had the same property since it was built. A run stops either way, naming a sample that
really failed.

## What was measured

All in the development container at this step's working tree.

| check | result |
|---|---|
| `cargo test --lib` | **5,905 passed, 14 ignored, 0 failed** — 5,885 at a2bf28c3 + 17 from the owner's merge + 3 from this step |
| `cargo test --lib ng::run` | **435 passed** (432 at the merged base + 3 added here) |
| `cargo test --lib ng::locus_generation` | 376 passed, 1 ignored |
| `cargo test --lib ref_seq` | 38 passed |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo doc --no-deps` | 26 unresolved-link errors, 23 redundant-link warnings — the standing baseline, unchanged |

**The three tests added:** `the_parallel_cover_gives_the_serial_drivers_answer` (the
handing-over driver against the undivided oracle, five widths × pools of 1, 2 and 8, on the
fixture carrying a boundary-straddling locus and a width-refused span),
`a_failing_source_ends_the_merge_under_the_parallel_cover`, and
`a_run_walker_can_cross_a_thread` (the compile-time `Send` proof). Beyond those, **every
existing record-path test now runs through the parallel cover** — the test harness's rayon pool
has more than one thread — so the fixtures that compare records, counts and the written VCF
against fixed expectations exercise the new schedule on every suite run.

**⚑ The base of this step moved mid-session, and the counts only add up once that is known.**
The session began at a2bf28c3, where `cargo test --lib ng::run` counts **421** (re-verified on
a clean checkout during this step). While E1 was being written, the owner merged the pre-pass
to calling handover branch into main (cce7d63d, 20 commits, 14 files — none of them files this
step touches). That merge adds 17 library tests, **11 of which match the `ng::run` substring
filter without being in `src/ng/run/`** — `ng::calling::run_parameters::…` contains the
substring `ng::run` inside `calling::run_parameters` — so the filter's count at the new base is
432, and 432 + the 3 tests this step adds = the 435 measured here. E1 is therefore committed on
top of cce7d63d, and the validation table below is against that base, not a2bf28c3.

## Deviations from the plan, recorded

1. **E1 does not build the pool of genotyping callers the step's text describes, and E2's
   "worker count" becomes the rayon thread count.** This is the milestone's own decision
   procedure — its deferral note says what E becomes is decided from the measurement, and the
   measurement is above. `CallersInFlight` stays a sketch in the architecture document, marked
   unbuilt with the number that ruled it out.
2. **Three types outside `src/ng/run/` changed**, all auto-trait plumbing with no API change,
   each at the site its own documentation had reserved (see *What landed*, point 3). The
   serial path's cost was measured before and after — see the checkpoint report's table — and
   the locks are uncontended by construction.
3. **19 `arc_with_non_send_sync` waivers were removed** (7 `#[expect]` that became compile
   errors the moment `WindowedRefSeq` became `Sync`, and 12 `#[allow]` that became dead text —
   counted from the patch: 19 deleted attribute lines, 7 under `#[expect(`, 12 under
   `#[allow(`; a first draft of this report said 17 and the artefact review counted it right).
   The one survivor, `examples/profile_posterior_engine.rs`, wraps production's
   `StreamingChromRefFetcher`, which is untouched.

## What the parallel cover bought, measured on the finished tree

Alternated arms in one sitting on a quiet machine, `NG_COVER=serial|parallel` on the probe,
63 samples, 8 rayon threads in the container:

| ground | serial calling | parallel calling | wall gain | serial CPU | parallel CPU |
|---|---|---|---|---|---|
| 200 kb (three pairs) | 12.06–12.32 s | 6.61–6.68 s | **1.8×** | 14.7–15.0 s | 22.2–22.8 s |
| the whole 8 Mb benchmark (one pair) | 473.6 s | 308.8 s | **1.5×** | 7m59s | 13m37s |

Identical loci on every run (23,450 and 1,069,772; 5,161 and 218,985 records). **The gain is
real and modest, and both gaps are understood rather than mysterious**: the parallel arm burns
about half again the CPU (Jacobi re-sweeps, per-cover scheduling, and records freed on a
different thread from the one that allocated them — the cross-thread-free cost Milestone G
exists to remove), and each cover is one 200-base building region whose sweep waits for the
slowest sample's decode, so the barrier fires some 40,000 times on the whole ground. Covering
a batch of regions per sweep, the way the merge's own parallel driver covers a whole round in
one call, is the next lever if the owner ever wants one — it changes only the cache's resident
window, not any answer. Recorded, not built.

## What the reviews found and what was done about it

Three agents reviewed the uncommitted step in isolated worktrees, each re-pointed at the base
commit with the step handed over as a patch: correctness with mutation testing, a determinism
reviewer whose only job was to break byte-identity across thread counts, and one whose only job
was to read what a person sees.

**The determinism attack found nothing, and what it threw is worth listing** because it is the
evidence behind "identical at every thread count": 432 parallel merges over 14-sample
deletion chains built to be the Jacobi worst case (a sweep per link), 240 over random layouts
at random pool sizes, the shipped fixture at 20 repetitions a pool, cache-level window
comparisons `cover_in_parallel` had never had directly, the two-failures race at 30 repetitions
a pool (always `Err`, always a sample that really failed — the accepted nondeterminism is which),
and real VCF files byte-identical across pools of 1, 2, 4 and 8. Its probe suite is saved under
the review directory and is the seed of E2's fixture.

**The correctness pass ran 14 mutations and killed 11; all 3 survivors are output-equivalent by
design**, not missed defects: the `> 1` pool switch flipped either way (both branches produce
identical output — the loss is scheduling, and nothing guards against silently falling back to
the serial cover, which is the same untestable-by-output property the parallel merge's own
module records), and the spare-recycling offer dropped (an offer, not an obligation, by the
trait's contract). The mutations that mattered were killed by the new test at its narrow widths
(`max→min` in the reduce), by a 150-second hang against a 7-second baseline (inverted fixpoint
break), and — the one worth quoting — dropping `covered_to` was killed by
`a_cohort_of_one_sample_writes_its_records`, which proves the one-sample cohort runs through
the parallel cover end to end. It also answered the deadlock question definitively: no path
re-enters a held lock, and the one drop that matters (the preparation guard before `shed`) is
compiler-enforced.

**The artefact pass found 7 Major and 8 Minor, all applied.** The two that mattered most: the
subcommand's `--help` and module doc still told a person the run is single-threaded — the one
behavioural change this step makes to the shipped command was the one the command denied — and
this report's own first draft carried two wrong claims about itself: the waiver count (17 where
the measured count is 19) and "three sites' documentation had reserved room" where only two
had. Both corrections are in the text above; the stale "5 ms per compressed megabyte" claim in
the probe's header, contradicted by its own new 63-sample row (9.84), is rewritten as the
measured range. Two stale `!Sync` mentions in owner-side documents
(`doc/devel/ng/spec/run_streaming.md` §8's trap entry among them) are recorded in
PROJECT_STATUS's overtaken-statements bullet rather than edited.

## What E1 does not do

- **No pool of genotyping workers, no `CallersInFlight` knob** — measured out, above.
- **`merge_cohort_in_parallel` is untouched and still unreached by a calling run** — it stays
  the merge's own oracle-checked parallel driver, off by default.
- **The walker is still not `Sync`**, so it still cannot go under the parallel merge's shared
  cache; `Send` is exactly what the sweep needs and no more.
- **Eviction stays serial in the calling drivers** — 0.8% of a run.
