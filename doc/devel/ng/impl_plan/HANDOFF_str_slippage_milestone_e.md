# Handoff — the STR slippage curve, Milestone E

*Written 2026-08-20 to start a fresh conversation. Paste the section below the line as the first
message; everything it needs is either in it or reachable from it.*

---

We are building **ng**, an experimental SNP/indel/STR caller; `PROJECT_STATUS.md` says where it
stands. Work on branch `ng-str-slippage-curve` in the worktree
`/Users/jose/devel/pop_var_caller-slippage-curve` — **not** the main checkout at
`/Users/jose/devel/pop_var_caller`, which another session is using. Always `cd` to the worktree
explicitly in every command; a drifting shell has already caused damage in this work once.

**Read `ai/skills/reporting-in-chat/SKILL.md` before your first reply and follow it in every
message to me** — it is not optional. `ai/skills/clear-technical-writing/SKILL.md` governs anything
you write into a document. `CLAUDE.md` at the repo root governs how you run things; two of its
sections bind everything you do: *What this caller has to work on — the range, not the example*,
and *Writing for the reader — including in chat*.

## The task

Execute **Milestone E** of
[`doc/devel/ng/impl_plan/str_slippage_level_curve.md`](str_slippage_level_curve.md), using the
`ai/skills/plan-driven-implementation` skill. The design is settled in
[`doc/devel/ng/spec/str_slippage_level_curve.md`](../spec/str_slippage_level_curve.md) §5.1 —
follow it rather than re-deciding it, and take any real design gap back to me rather than closing
it yourself.

**In one sentence:** the slippage *level* already gets a curve across repeat count with each
stratum departing from it by how much evidence it has; Milestone E gives the **direction split**
and the **fall-off** the same treatment, and retires the gate-and-copy rule that currently serves
them.

## What is already built and committed

Milestones A–D are done. In `src/ng/parameter_estimation/joint/`:

- **`slippage_curve.rs`** — the curve for the level. `RiseShape` (0 = each repeat multiplies the
  level, 1 = adds), `fit_line`, `choose_rise_shape` (grid search scored by leaving each cell out),
  `blend_level` (inverse-variance on the log scale, with a knee that stands the curve down where
  it disagrees with a cell by more than either error explains), and the provenance types.
- **`ssr_fit.rs`** — `fit_strata` now fits every stratum on its own tracts, copies the two shares
  from a neighbour where its own slipped reads fall short of 4,000, draws a curve a motif period,
  and re-emits every level through the blend. **Nothing pools tracts any more.**
- **`census.rs`** — `RECORDED_OFFSET_RANGE` moved 4 → 8 on measurement (below).

Tests: 36 in `slippage_curve.rs`, 11 in `ssr_fit.rs`, 5 in
`tests/slippage_curve_on_real_cells.rs` which run the fit over both cohorts' real cell tables under
`tests/data/slippage_cells/`.

## What Milestone E must do, and why

**D's rule delivers almost nothing, and that is the measurement that motivated E.** With the
floor-and-copy rule as built:

| | strata fitted | furnished from a neighbour | still refused |
|---|---:|---:|---:|
| HG002, one sample at ~300 reads | 55 | 13 | 69 |
| tomato, 63 accessions at ~3 reads | 6 | **0** | 65 |

because the best-measured stratum of each motif period holds this many slipped reads, against the
4,000 a stratum needs before its shares count as its own:

| best slipped reads at | period 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| HG002 | **8,840** | 2,562 | 314 | 212 | 39 | 0 |
| tomato | 2,338 | 128 | 0 | 0 | 0 | 0 |

Only HG002's homopolymers clear it. Spec §5.1 replaces the gate with a curve fitted from *every*
stratum weighted by its precision.

**Three things the shares need that the level did not**, all in spec §5.1:

1. **A family that can bend either way, and it differs by period.** The level's family is
   monotone and refuses a falling fit; the shares do neither. Measured at ±8 from strata that
   copied nothing (report §5): HG002's dinucleotide direction split spans **4.52-fold** and a
   logit-line beats the mean four to one (0.060 against 0.240); its dinucleotide fall-off spans
   5.17-fold and **nothing beats the mean**; the homopolymer split rises 0.51 → 0.82 while the
   dinucleotide one falls then climbs. E1 chooses among a constant, logit-linear and
   logit-quadratic **by the held-out-cell criterion**, per period and per parameter. *A period
   whose held-out error is lowest at the constant has no trend to fit, and that is also the answer
   to the research plan's C4.*
2. **A share's own precision:** `sqrt((1 − p) / (p · S))` on `S` slipped reads. This reproduces the
   architecture's own two numbers — 1,357 slipped reads to hold a split of 0.17 to 6% where it
   says "about 1,400", 3,997 for a fall-off of 0.065 where it says "about 4,000" — so it is the
   model the 4,000 came from.
3. **Thin strata must be fitted to contribute at all.** Nothing below 50 tracts is fitted today.
   E4 lowers that, **gated on** spec §11's open question: fit drawn strata down to a handful of
   tracts first and report how often the climb converges and how often the level comes back
   exactly zero. If those fits are unusable, E4 does not happen and the floor stays.

**One rule that must not be dropped:** a stratum feeds its curve through its **own** fit, never a
blended one, or each round of smoothing fits a curve to the previous round's curve.

## Two traps that have already cost time here

**The census's recording window was lying, and the report has been corrected for it.** It
recorded a read's length offset over only ±4 repeats; at 30-repeat homopolymers that under-measured
the level by **2.26×**, and it made the rise look as though it flattened above 20 repeats when it
does not. The window is now ±8, converged against a ±12 arm agreeing within 1.8%.
[`../reports/str_slippage_shape_2026-08-20.md`](../reports/str_slippage_shape_2026-08-20.md) has
been remeasured throughout, and its **§7 lists what the correction overturned — read that before
trusting any older number you find elsewhere in the docs.** The ±8 cell tables are in the worktree
at `tmp/slippage_curve/hg002_d4*.csv`, `hg002_wide8_plain.csv` and `tomato_d4*.csv`. **For the two
shares use `hg002_wide8_plain.csv`, not `hg002_d4_plain.csv`** — the second has shares copied from
neighbours in it, and a copied share read as a measurement is exactly the circularity Milestone E
exists to avoid.

**`SSR_BORROWING_FLOOR` no longer exists** and the env knob is now `SSR_SHARES_FLOOR`. An earlier
run of mine set the old variable, got zero furnished strata, and I reported that as a design gap
when it was a run setting. Check what a switch actually controls before concluding from it.

## How to run things

- Builds and tests: `./scripts/dev.sh cargo …` from inside the worktree (Apple `container` on this
  machine). **Never `cargo fmt` repo-wide** — format the files you touched:
  `./scripts/dev.sh cargo fmt -- <path>`.
- The aggregate `cargo clippy --all-targets` gate is red and pre-dates this branch, in
  `src/ng/run/cohort_merge/observation_cache.rs` and `examples/shared/synthetic_alignment.rs`.
  `cargo clippy --lib --all-features -- -D warnings` is clean and is the gate to hold.
- A cohort run needs the main checkout mounted read-only for its data:
  `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller ./scripts/dev.sh env … ./target-container/release/examples/ng_joint_records_walk <ref> <catalog.parquet> <regions.bed> <generic-target> <alignments…>`
  with `SSR_CELL_TABLE=` and `SSR_CELL_TABLE_NO_CURVE=` for the two arms. HG002 takes about 20
  minutes, tomato about an hour. The exact invocations are in
  `tmp/slippage_curve/*.log`'s first lines.
- Scratch goes in the worktree's `tmp/`, never the system `/tmp` and never the harness's own
  scratch directory.

## What I want from you

Work step by step and **pause at each milestone checkpoint**. When you report, give me a
recommendation rather than a list of options, and never assert a property without its size, its
subject and its measure. E1 produces a report and no code; bring me its numbers before E2 starts.
