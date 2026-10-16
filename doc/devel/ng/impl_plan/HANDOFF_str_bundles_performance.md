# Handoff: performance of the STR bundle caller

*2026-10-15. For a fresh session that acts on the performance review. Read this, then the review:
[`../../reports/reviews/str_bundles_performance_2026-10-15.md`](../../reports/reviews/str_bundles_performance_2026-10-15.md).*

---

## 1. Where things are

- **Worktree** `/Users/jose/devel/pop_var_caller-bundle-perf`, **branch** `bundle-perf`, from main e4bd1e94.
  Nothing in `src/` is changed yet; the branch holds the review and this note.
- **The feature under review** is the STR bundle caller: two or more repeats within 15 bases of each
  other called as one STR locus, with stutter in each pure run
  ([`../spec/str_bundle_loci.md`](../spec/str_bundle_loci.md); results in
  [`../../reports/implementations/str_bundles_2026-10-12.md`](../../reports/implementations/str_bundles_2026-10-12.md)).
- **The container service** (Apple `container`) was down on 2026-10-15; the review built and measured
  **natively** (`cargo build --release`, binaries in `target/release`). That is fine for measuring and
  profiling — macOS profilers (`samply`) work natively — and the owner allows native cargo in the dev
  loop. The full test suite still goes through `scripts/dev.sh` once the owner restarts the service
  (`container system start`; auto mode refuses to run it). In the container use `cargo test -j 2`: a
  16 GB container was OOM-killed building the test binary with the default parallelism.

## 2. The owner's brief

> Targets: wall time and memory. Measure rather than infer. Act on the review and check whether each
> change works. A change that makes the code clearer, simpler or better structured is accepted even with
> a small gain or a small loss. A change that complicates the code for a small or no gain is rejected.
> Clear trade-offs (half the wall time for 1% more memory) are wins; hard ones are discussed with the
> owner.

## 3. The baseline (native, M5 Pro, 18 cores)

Bundle cost is the difference between the default and `--max-str-bundle-len 0` (every bundle a
satellite, uncalled).

| workload | command | threads | wall ON / OFF | instructions ON − OFF | peak RSS ON / OFF |
|---|---|---|---|---|---|
| 50 samples, 30× | generate-psps | 1 | 322 / 291 s | +558 G of 9,199 G (6.5%) | 189 / 190 MB |
| 50 samples, 30× | call-from-psps | 1 | 36.7 / 33.7 s | +75 G (13%) | 262 / 261 MB |
| 50 samples, 30× | call-from-psps | 8 | 8.6 / 7.2 s | +75 G | 288 / 284 MB |
| 50 samples, 5× | call-from-psps | 1 | 17.0 / 16.2 s | +24 G | 261 / 258 MB |
| HG002 30× | call-from-alignments | 1 | 33.1 / 30.5 s | +54 G | 287 / 287 MB |
| HG002 30× | call-from-alignments | 8 | 30.7 / 28.2 s | +53 G | 285 / 280 MB |

Memory: no finding (ON and OFF within 11 MB). Cost is linear in samples.

## 4. The findings, in the order to try them

All must keep the output byte-identical (VCF bodies and psps). Clarity effect in brackets.

1. **`render_tract`: size the buffer once, keep renderings in one flat buffer** — up to 18.5% of bundle
   calling (27% at 5×); helps segments too. [simpler]
2. **Reuse the first run parse for admitted alleles** instead of parsing twice. [simpler]
3. **`parse_into_runs`: store the table in fill order, flatten the template, reuse scratch** — 23% of
   bundle calling. [neutral]
4. **One STR read cursor for segments and bundles** — about 77% of bundle generation (≈25 s CPU over 50
   samples at 30×; ≈1.5 s HG002 wall at any thread count). The largest gain, at a moderate cost in
   structure: measure it and weigh it against the owner's rule; this is the one most likely to need the
   owner's call.
5. **Lay each candidate out once for complete reads** (26% of bundle calling at 30×): buffer reuse first
   [neutral]; the full cache only if the measured gain is large [more complex].
6. **Cache partial-read layouts per allele, probabilities per read group** (14–27%). [more complex]
7. Smaller: two logarithms per rendering; walk inner-loop tweaks and slip-set bounds; index jumps past
   long gaps (cheap alternative to 4); noodles-bgzf `libdeflate`; 2–4 merge regions per worker (outside
   bundle code).

**Outside the bundle scope, noted by the review:** `memmove` is 19% of `call-from-psps` self time, in the
psp read path — more than all of bundle calling. It deserves its own review; raise it with the owner.

## 5. How to measure each change

Scripts in `/Users/jose/devel/pop_var_caller-bundle-perf/tmp/perf/` (not committed; they live in the
worktree's `tmp/`):

- `bundle_baseline.sh [sim] [hg002] [table]` re-runs the baseline and prints the table; `TAG=<name>`
  gives a variant build its own results file. Baseline results: `results_base.tsv`.
- `profile_bundles.sh` re-takes the profiles; `samply_rank.py` ranks a profile headlessly;
  `profile_diff.py` subtracts an OFF profile from an ON one; `table.py` prints the table.
- Gate each change on **instructions retired** (stable on a loaded host; wall time is not) and
  **byte-identical output**; report wall and RSS too. Record `uptime` with every measurement — other
  sessions share this machine and have pushed its load past 100.
- Data (read-only): `/Users/jose/devel/pop_var_caller-str-bundles/tmp/bundle_bench/runs/{pcr_50s_30x,
  pcr_50s_5x}` (BAMs, reference, catalog, true stutter); HG002 under
  `/Users/jose/devel/pop_var_caller/benchmarks/`.
- `tmp/perf/scratch/` holds 5.6 GB of baseline psps; clear it once no change needs them.

## 6. Workflow the owner expects

One agent implements each change with `ai/skills/rust-feature-implementation/SKILL.md`, measuring
before/after; a second agent reviews and fixes with `ai/skills/rust-code-review/SKILL.md` and
`ai/skills/apply-code-review-fixes/SKILL.md` (no per-category git worktrees — disk). Report to the owner
with `ai/skills/reporting-in-chat/SKILL.md`: answer first, numbers with their denominators, one
recommendation per decision. Commit on the branch; merging into main needs the owner's yes.

## 7. Open items from the bundle work, not about performance

- **Main has not been through the full test suite since the merges of 2026-10-13** (simulator hard
  shapes, 206a4832; partial reads nominating candidates, 0e905bbb): the container service was down. Run
  it once the service is up.
- **The owner was offered a whole-feature code review** of the bundle caller (each step was reviewed
  alone; the merge e7c9106e that reconciled the segment fix with bundles never was). Not yet answered.
- **Dropped, kept unmerged:** branch `bundles-norun-generic` (bundles with no run → generic path).
