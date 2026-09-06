# ng — window coverage: implementation plan

*Draft, 2026-09-06. Turns [`../spec/window_coverage.md`](../spec/window_coverage.md) into build
order — no new design here; a design question surfacing mid-plan goes back to that spec. There is
no architecture document: the spec's §3 type blocks are the code shape the steps cite. The
consumer of what this builds is [`hidden_paralog_filter.md`](hidden_paralog_filter.md), which
starts after this plan's Checkpoint C. **The standing oracle for every milestone is that a run
with the filter off writes byte for byte what it writes today** — on the run fixtures and on the
tomato slice ([`scripts/ng_mode_equivalence_oracle.sh`](../../../../scripts/ng_mode_equivalence_oracle.sh)
produces both modes' files). A step that moves a byte is wrong, not almost done.*

## Scope

**In:** the sliding-window accumulator transcribed from production with its tests, plus the floor
and the per-sample depth scale; the per-position depth rule over `Drawn` (head at single-base
records, body elsewhere, anchor-only for generic records); the reference bases into the cache;
the accumulator in the cache's per-sample state, observed at the draw, held and evicted with the
summaries; the half-window look-ahead on `cover`; the pair on the cohort locus and on the output
evidence, and beside the record at the sink; the histograms out of the cache; the whole-store
recomputation probe; the three measurements the spec leaves to this plan.

**Out (later plans):**

- **anything that reads the pair or the histogram** — the model fit, the score, the spill, the
  verdict: [`hidden_paralog_filter.md`](hidden_paralog_filter.md);
- **a histogram in the psp trailer** — spec §8, only if a measurement forces it;
- **depth for `SsrBundle` and `Satellite` ground** — [`locus_generation.md`](../spec/locus_generation.md) §11;
- **releasing the psp source's arena** — [`cohort_merge_psp_path.md`](cohort_merge_psp_path.md)'s
  own first debt; this plan neither fixes nor worsens it.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone.
- **The copy before the change.** Production's accumulator is transcribed *verbatim* with its
  tests and proven green before the floor or the depth scale touch it — the same
  transcribe-first rule the pileup port paid for three times
  ([`locus_generation/pileup/mod.rs`](../../../../src/ng/locus_generation/pileup/mod.rs)).
- **The algorithmic heart before the plumbing.** The depth rule over a `Drawn` is a pure
  function with its own tests before the cache learns to call it.
- **Verify against ground truth, not self-consistency.** The oracle is a probe that recomputes
  every window from a whole store independently of the cache (spec §10); the cache's values must
  equal it bit for bit.
- **Isolate the silent step.** The look-ahead on `cover` fails by leaving every region's last
  centres absent — no panic, no test failure unless one looks — so it lands alone, oracle green
  before and after.
- **Measure, don't guess.** The floor, the bin scheme and the memory are numbers this plan
  produces; the spec marks its starting values soft for that reason.
- **Container builds**: `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (verify before step A1)

- **Both modes reach the merge through one function** —
  [`call_cohort_from_sources_handing_each_record_over`, `callers.rs:799`](../../../../src/ng/run/callers.rs),
  with psp mode building `PspSummarySource` ([`psp_caller.rs:477`](../../../../src/ng/run/psp_caller.rs))
  and the cache materialising a kept body at assembly through `source.build`
  ([`observation_cache.rs:710`](../../../../src/ng/run/cohort_merge/observation_cache.rs)). On
  `main` as of 2026-09-06; confirm by running the mode-equivalence oracle green before A1.
- **The psp record head carries `reads_compared_with_reference`** —
  [`psp_head_compared_reads.md`](psp_head_compared_reads.md) Milestone H, met.
- **`num_obs_along_locus`** ([`locus_generation/mod.rs:67`](../../../../src/ng/locus_generation/mod.rs))
  spreads partial witnesses per run — the C2 change its comment describes is in.
- **Production's accumulator and its tests** are where the spec says
  ([`coverage.rs:331-1050`](../../../../src/sample_summary/coverage.rs)), and production is
  frozen: nothing in `src/sample_summary/` is edited.
- **A psp store to measure on**: the six tomato accessions over `tmp/c1_two_regions.bed`
  written by [`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh),
  and one whole-genome tomato psp (`SRR7279481`, the store the skip measurement was re-taken on
  2026-09-04) for the single-base equality at scale.

---

## The steps

### Milestone A — the accumulator, copied and then extended

**A1. ✅ Transcribe the accumulator.** `src/ng/window_coverage/mod.rs` (the config, `WindowCoverage`,
the histogram type) and `accumulator.rs`: `SlidingWindowCoverageAccumulator`, `CoverageBinScheme`,
`WindowCoverage`, `CoverageByGcHistogram`, and the sliding-window tests (`sliding_*`,
[`coverage.rs:780-1050`](../../../../src/sample_summary/coverage.rs)) transcribed unchanged;
the fixed-tile accumulator and the heterozygosity fields are not brought over. Green as
transcribed. *Depends:* —. *Source:* spec §3.3, §3.4, §7.

**A2. ✅ The floor.** `min_window_positions` on the config; a centre finalised over fewer covered
positions emits `NaN` in both fields. Tests: a window at the floor emits a value, one below it
emits `NaN`, and the `NaN` compares by bits. *Depends:* A1. *Source:* spec §3.3.

**A3. ✅ The per-sample depth scale. Own commit, do not bundle.** The first
`depth_scale_windows` finalised windows are buffered; at the 10,000th, the width is set from
their median and every buffered window is folded; `finish` handles fewer than 10,000 and none.
Its failure is silent — a width off by a factor is a histogram whose mode sits in the wrong bin,
and every fit anchors on it — so the tests pin: the width from a known stream; the same
histogram from a stream folded early and late; a stream shorter than the scale sample; an
empty one. *Depends:* A1. *Source:* spec §3.4.

> **Checkpoint A: production's window, green under production's tests, with ng's two
> additions pinned. Pause for review.**

### Milestone B — the depth at a position

**B1. ☐ The depth rule over a `Drawn`.** A pure function: given a `Drawn` and, for a kept
record spanning more than one base, the source's `build`, yield the positions the record speaks
for and the depth at each — the head count at a single-base record; the body's
`num_obs_along_locus()` at the anchor only for a generic record, at every position for a tract.
Tests on fixtures of each shape, both `Built` and `Kept`. *Depends:* —. *Source:* spec §3.1.

**B2. ☐ The single-base equality, measured. Own commit.** A probe (`examples/ng_window_coverage_probe.rs`)
walks a real store building every body and asserts, at every record spanning one base, that the
head count equals the body's `num_obs` sum; reports the count checked and any counter-example,
on the six-accession slice and the whole-genome `SRR7279481` store. **A counter-example changes
spec §3.1 to "build every body" before C begins.** *Depends:* B1. *Source:* spec §3.1, §6 trap 7.

> **Checkpoint B: the rule is a function with tests, and the head-equals-body claim is a
> number on a real store. Pause for review.**

### Milestone C — into the cache, and out to the record

**C1. ☐ The reference into the cache.** `ObservationCache::over` takes a reference accessor
minted beside the run's padding one (`walk_reference.accessor()`, both callers); `cover` fetches
the cover's ground once into a scratch buffer with `fetch_into`; a fetch failure is a `RunError`
naming the region. Direct mode's VCF byte-identical. *Depends:* —. *Source:* spec §3.2.

**C2. ☐ The accumulator in the per-sample window.** `SampleWindow` gains one accumulator and a
deque of finalised `(position, WindowCoverage)`; `draw_to` observes each drawn record through
B1's rule against the fetched bases and drains `pop_ready` into the deque; `evict_before` and its
parallel twin drain the deque with the summaries; `WindowedCohort` gains a third view over it,
`window_at(sample, position) -> Option<WindowCoverage>`. The parallel cover's per-sample `&mut`
covers the new state. VCF byte-identical. *Depends:* A3, B1, C1. *Source:* spec §3.3, §3.5.

**C3. ☐ The look-ahead. Own commit, do not bundle.** `cover` draws to the region's reach plus
half a window. The probe from B2 is extended to record every position's pair by spec §3.3's rule
independently of the cache (the whole-store recomputation, spec §10); a test-only hook reads the
cache's deque at every built locus and asserts bit-identity with the probe, on the six-accession
slice. **Green before and after this commit**, so that the absent-at-region-edges failure is
visible if it ever returns. *Depends:* C2. *Source:* spec §3.3, §6 trap 3.

**C4. ☐ The pair on the locus and beside the record.** `CohortObservation` gains
`window_coverage`, parallel to `per_sample`, read at the locus's first base;
`SampleEvidenceForOutput` gains the dense pair; the sink both callers take becomes
`FnMut(&VcfRecord, &[WindowCoverage]) -> Result<(), E>`, with the VCF writer's adaptor
ignoring the slice. VCF byte-identical; the mode-equivalence oracle extended to compare the
pairs, green. *Depends:* C2. *Source:* spec §3.5.

**C5. ☐ The histograms out.** `ObservationCache::into_sources` is joined by a form that finishes
every accumulator and returns the histograms beside the sources; both callers carry them to
their outcome type, unused until the filter plan. The probe's histogram equals the run's, byte
for byte, per sample. *Depends:* C2. *Source:* spec §3.4, §3.5.

> **Checkpoint C: every written record carries the pairs, every sample ends the pass with its
> histogram, both equal the whole-store recomputation bit for bit, both modes agree, and no
> VCF byte moved. The filter plan may start. Pause for review.**

### Milestone D — the three numbers the spec left to measurement

**D1. ☐ The floor.** From the probe: the distribution of covered positions per window on the
tomato slice and on HG002; what share of windows fall under 50, 100, 200; the default written
into the config with the distribution beside it in the report, and spec §3.3's OPEN closed.
*Depends:* C3. *Source:* spec §3.3, §9.

**D2. ☐ The bin scheme.** From the probe: each sample's fitted width, and the overflow fraction
— windows past the range — on both benchmarks, against the fit's rejection guard at a fifth; the
scale sample and the factor of ten confirmed or moved, spec §3.4's OPEN closed. *Depends:* C5.
*Source:* spec §3.4, §9.

**D3. ☐ The memory, per sample.** Peak resident of a psp-mode run at 1, 6 and 63 samples on the
tomato slice, before and after this plan, with the per-sample slope reported against
[`run_streaming.md`](../spec/run_streaming.md) §7.2's budget; the histogram's 80 kB and the
look-ahead's held summaries priced separately. *Depends:* C5. *Source:* spec §3.4, §4, §5.

> **Checkpoint D: the floor and the bins are measured defaults, and the cost is a number
> against the budget. Pause for review.**

## Verification summary

| milestone | proven by |
|---|---|
| A — the accumulator | production's sliding-window tests green as transcribed; the floor and the scale pinned by their own |
| B — the depth rule | fixture tests over both `Drawn` shapes and all three record kinds; the single-base equality asserted over a real store |
| C — the cache | VCF byte-identical at every step; the pairs and histograms bit-identical to the whole-store recomputation; mode equivalence with pairs compared; thread-count invariance |
| D — the numbers | the three measurements reported with their distributions, on both benchmarks |

## Out of scope (next plans)

- **[`hidden_paralog_filter.md`](hidden_paralog_filter.md)** — everything downstream of the
  pair and the histogram.
- **The trailer histogram** — [`../spec/window_coverage.md`](../spec/window_coverage.md) §8; no
  plan until a measurement asks for one.
- **Bundle and satellite depth** — [`locus_generation.md`](../spec/locus_generation.md) §11's
  generator.
