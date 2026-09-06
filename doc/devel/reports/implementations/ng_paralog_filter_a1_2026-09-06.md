# The hidden-duplication filter — A1: the constants and the coverage model, copied into ng

**Date:** 2026-09-06
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone A, step A1
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.1, §7
**Branch:** `ng-paralog-filter`

## The answer

**Production's per-sample coverage model now lives in `src/ng/paralog/`, byte for byte, and a
test says so.** The model answers one question — *at a window of this GC content, what read depth
would one copy of this sequence produce in this sample?* — and the filter divides an observed
window depth by its answer to get a relative copy number. Two copies of a gene the reference
collapsed into one read as about two.

The copy is not merely tested, it is **diffed**: [`copy_fidelity.rs`](../../../src/ng/paralog/copy_fidelity.rs)
compares each of ng's two copied files to its production original — line for line and then byte
for byte — and fails the build on any difference. Production's own tests, transcribed with the
files, would stay green through a divergence — both trees have their own copy of them — so a
green suite is not evidence of fidelity and a textual comparison is.

**The review this step went through found that the guard covered less than its prose said, in
three ways that each let a real divergence through**, and all three are closed and proved closed
by re-running the mutation that used to pass: a drift in the copied constants (they were in the
one file the guard exempted), a guarded-file list that could go empty while still printing `ok`,
and a stripped final newline or a `CRLF` ending, which a line comparison cannot see. See the
[review](../reviews/ng_paralog_filter_a1_2026-09-06.md) and the
[fixes applied](../reviews/fixes_applied_ng_paralog_filter_a1_2026-09-06.md).

## What was added

| file | what it is |
|---|---|
| [`src/ng/paralog/coverage_model.rs`](../../../src/ng/paralog/coverage_model.rs) | production's file, byte for byte, with a note appended to its module header |
| [`src/ng/paralog/model_params.rs`](../../../src/ng/paralog/model_params.rs) | production's `src/paralog/mod.rs` from its first item to its end — the constants and grids (`ParalogModelParams`, `SfsPriorSpec`, `GridSpec`) and the four tests that pin them — under an ng header |
| [`src/ng/paralog/mod.rs`](../../../src/ng/paralog/mod.rs) | ng's own: the module declarations, the re-exports, and one test of ng's beside production's |
| [`src/ng/paralog/copy_fidelity.rs`](../../../src/ng/paralog/copy_fidelity.rs) | ng's own: the textual guard, its own fixtures, and the record of which files are still copies |

`src/ng/mod.rs` gains `pub mod paralog;` and a clause in its module header. Nothing under
`src/paralog/` or `src/var_calling/` was touched.

**The constants live in their own file so that the guard can reach them.** They were first put
in `mod.rs` beside ng's declarations, where a file-level comparison cannot go — production keeps
them in a `mod.rs` whose declarations ng does not share. The review proved the cost: deleting
both `is_finite()` guards from `GridSpec::new` passed all 6,303 lib tests, `copy_fidelity`
included. Production keeps them in a `mod.rs`, so this copy is compared from the first item's
`///` on each side rather than from the end of a shared module header.

**The copies were diffed, not assumed.** `coverage_model.rs` matches production's over
**1,158 lines** past the header, and `model_params.rs` matches `src/paralog/mod.rs`'s items over
**236 lines**. Both diffs are empty.

The constants are the tomato2 prototype's and are **inherited, not re-measured** — the
single-copy band `0.4..1.6`, the minimum bin count 50, the five-bin smoothing window, the
mode/median guard `[0.5, 1.5]`, the overflow guard at a fifth (spec §3.1), and on the model side
the carrier copy numbers `{3, 4, 6, 8}`, the winsor cap at four copies, the 200-point allele
frequency grid, the 40-point carrier frequency grid on `[0.004, 0.6]`, the `0.01` error floor and
the hom-alt veto at VAF 0.9 with five reads (spec §3.2).

## Assumptions and recorded deviations

**1. The histogram the fit reads is production's type, and it is meant to be ng's.** The spec's
reuse map (§7) says the copied fit "reads ng's histogram type, which is production's shape", and
window coverage's spec §3.5 puts ng's `CoverageByGcHistogram` in `src/ng/window_coverage/`.
**That module is not on this branch or on `main`**: it is on `ng-window-coverage`, which had no
commits when this step began. So the plan's first precondition does not hold — and it fails one
step earlier than the plan says, because this step needs window coverage's *first* step, not its
Checkpoint C.

The copied fit therefore still names `crate::sample_summary::CoverageByGcHistogram`, exactly as
production's does. The *fit* reads five of that struct's fields — `window_bp` once, `gc_bins`
twice, `depth_bin_width` six times, `depth_bins` once, `counts` once — and nothing else.

**The swap will cost four lines, not one, and the review is what found that.** Checked against
the other branch's code rather than its spec: `ng-window-coverage` committed `09bad846` ("wip: A1
for review") while this step was under review, and its `CoverageByGcHistogram` carries
`window_bp`, `gc_bins`, `depth_bin_width`, `depth_bins`, `windows_folded` and `counts` — all five
the fit reads. But **the file does not end at the fit**: its transcribed test fixture
(`coverage_model.rs:714-723`) builds the histogram with all eight of production's fields,
including `n_positions`, `n_skipped_tiles` and `callable_positions`. So the swap is the `use`
line plus three lines of that fixture, inside a file the guard forbids editing — which is why it
has to land as its own commit that releases the file and says what it did.

**What that leaves open is another plan's decision, raised at Checkpoint A, not taken here.**
Window coverage's spec §7 says the fields the model fit does not read are dropped; production's
transcribed tests build all eight. Both cannot hold. Which way it resolves belongs to that
document.

**2. The fidelity guard is an addition the plan does not ask for.** The plan asks for the files
"green as transcribed". A textual guard is the codebase's existing answer to the same problem
([`locus_generation/pileup/copy_fidelity.rs`](../../../src/ng/locus_generation/pileup/copy_fidelity.rs)),
and it is the only check that can tell a transcription slip from a deliberate change. Its rule
here is narrower than the pileup port's in one way — these files need no path repoints, so the
sanctioned-exception machinery is dropped rather than carried empty — and wider in three: the
comparison is byte-level as well as line-level, it has its own fixtures rather than only the real
file pair, and the module's directory listing is checked against the guarded set so the guard
cannot go quietly empty. All three came out of the review.

**3. `mod.rs` holds nothing ng did not write, and that took a second file.** The constants were
first put in `mod.rs` beside ng's declarations, on the reasoning that its four transcribed tests
pin the values. The review showed that is not enough: those tests pin *values*, not code, so
deleting both `is_finite()` guards from `GridSpec::new` left every one of them green. The block
moved to `model_params.rs`, where the guard reaches it.

## Tests

31 tests under `ng::paralog::`, of which 27 are production's transcribed and four are ng's own.

**Production's, transcribed — 23 in `coverage_model`**, of which **seven are rejections**:
`rejects_empty_histogram`, `rejects_bottom_depth_bin_only_sample`,
`rejects_wrong_copy_number_peak`, `rejects_when_depth_range_exceeded`,
`rejects_when_all_gc_bins_thin`, `rejects_invalid_config` and `rejects_all_invalid_config_shapes`.
Seven more are the fit itself: its happy path with `relative_copy_number`, invariance to a scaled
count matrix, the GC curve flattening a planted GC shift, σ₀ from a known spread, σ₀'s floor when
the single-copy band collapses, the σ₀ override, and a small overflow tail being accepted rather
than rejected. The remaining nine are the helpers — two weighted-median cases, two gap-fill
cases, the smoothing, two mode-refinement cases, the GC multiplier's clamp beyond the bin
centres, and the mode/median bound constructor.

**Production's, transcribed — 4 in `model_params`**: the default parameters against the prototype
constants, the two grid defaults, the relation that the top carrier's coverage mean equals the
winsor cap, and `GridSpec::new`'s rejection of every degenerate grid.

**ng's own — 4.** Three in `copy_fidelity`: the copies are still production's; the comparison
itself accepts an appended module-header note and rejects every other edit, driven on synthetic
file pairs because the real pair is identical and would otherwise exercise each rejecting branch
in its accepting direction alone; and every `.rs` file in the module is guarded, ng's own, or
explicitly released. One in `mod`: `grid_spec_new_refuses_an_infinite_endpoint`, because the
transcribed `// non-finite` case does not reach the finiteness guard — `NAN < 0.6` is `false`, so
that input is refused by the `lo < hi` clause with or without `is_finite`, and an *infinite*
endpoint is what separates them.

**Each of the guard's three additions was proved to fire** by re-running, in this tree, a
mutation that passed before it: a drifted doc line in the copied constants, a stray `.rs` file in
the module, an emptied guarded set, and a stripped final newline. Every mutation was restored
from a backup and the restoration confirmed by `diff`.

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::paralog::"
    → test result: ok. 31 passed; 0 failed; 0 ignored; 6290 filtered out

    cargo test --all-features --lib --bins --tests
    → 6306 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

**That failure is `main`'s, and this step's diff was removed to prove it.** With
`pub mod paralog;` taken back out of `src/ng/mod.rs` — so that none of this step's files are
compiled at all — the same command gives **6275 lib tests passed** and the *same* single
integration failure, `tests/ng_calling_loop_calls_genotypes.rs:1235`,
`left: [0, 0]` against `right: [0, 1]`. So this step adds 31 tests and moves nothing:
6275 → 6306.

## What is red on main, and how this step was gated instead

Four checks fail on `main` at `a33ada0f`, on files this step does not touch:

1. **`cargo test --all-targets --all-features` does not compile.**
   `examples/ng_candidate_selection_probe.rs` reads `locus.kind` and hands a `ClosedLocusRanges`
   to `CohortObservation::over`, which takes a `ClosedLocus`. Three errors; no test runs at all,
   so the command exits 101 having measured nothing. The type it uses arrived at `d44e4eaf` on
   another line of work, after the example was last touched at `c90169ce`.
2. **One integration test fails** — the contamination-at-a-tract case above.
3. **`cargo fmt --check` is dirty on nine files** — `examples/ng_call_cohort_end_to_end.rs`,
   `examples/ng_call_from_psps_cost.rs`, `src/ng/psp/block.rs`, `src/ng/run/psp_source.rs` and
   four under `src/ng/run/cohort_merge/`.
4. **`cargo clippy --lib --bins --tests -- -D warnings` fails on three lints**, all
   `needless_lifetimes`: `cohort_merge/build.rs:820`, `:893`, `cohort_merge/serial.rs:67`.

None is this step's to fix, and fixing them would put twelve files belonging to other work into
this branch's commits. The gate used instead, and for every step of this plan until `main` is
clean: **`cargo test --all-features --lib --bins --tests`, with no failure `main` does not
already have**, and the `fmt` and `clippy` failure sets no larger than `main`'s — nine files and
three lints, none of them ng's paralog files.

## Follow-ups

- **The histogram swap when window coverage lands** — the `use` line in `coverage_model.rs` plus
  three lines of its transcribed fixture, released from the fidelity guard in its own commit.
  Which fields ng's histogram carries is window coverage's spec's call; raised at Checkpoint A.
- **Nothing yet makes the temporary import temporary** (review Mi3). The fix the review suggests
  — a step in the plan with window coverage's first step as its precondition — is a plan edit
  this loop may not make; recorded in `PROJECT_STATUS.md` and raised at Checkpoint A.
- **`main`'s four red checks** — not this plan's, but they cost every branch its `--all-targets`
  gate, and one of them is a genotype changing at a repeat tract.
- **`grid_specs_have_expected_defaults` reads as duplication** of its neighbour except for one
  assertion (it reaches `SfsPriorSpec::default()` directly). A sentence saying so belongs on the
  test, which is production's and guarded, so it is recorded here instead.
