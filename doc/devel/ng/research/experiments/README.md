# Archived experiments

Working documents from measurement work whose conclusions are written up elsewhere. They are
kept because they carry two things a conclusion cannot: the numbers that were rejected, and
code that was built, measured and then deliberately not adopted.

**The code in the diffs here is throwaway.** It was written to be measured, not to be merged.
Do not treat any of it as a starting point without re-reading why it was set aside.

---

## `locus_stream_shape/` — how the calling pipeline should pass loci between stages

Four experiments, run 2026-08-06, that settled whether a locus travels through the pipeline as
one self-contained object or many loci travel as one array per field. **The conclusion is
[../locus_stream_shape.md](../locus_stream_shape.md)**; the plan they were run under is
[../../impl_plan/locus_stream_shape_experiments.md](../../impl_plan/locus_stream_shape_experiments.md).

Each experiment kept the numbering it was given while running; this table is what each one
actually measured.

| files | what it measured |
|---|---|
| `sketch1_generator_output.md`, `sketch1_producer.diff`, `sketch1_src_block.rs`, `sketch1_example_ng_prepass_sketch.rs` | The pileup filling blocks of arrays instead of emitting one object per covered base, with the parameter pre-pass reading it three ways. **Worth 7.8 % of a covered base at 30× for 830 lines in the walk's hottest files** — the deferred decision. |
| `sketch2_merge.md`, `sketch2_code.diff` | The k-way merge across samples, reading a `.psp` file and reading the pileup directly, with and without a first pass over a cheap summary column. |
| `sketch3_calling_input.md`, `sketch3.diff` | Whether the calling step cares what shape its input is. It does not: the arithmetic is 167× the data handling. |
| `sketch4_columnar_producer_plus_fold.md`, `sketch4_code.diff`, `sketch4_src_block.rs`, `sketch4_example_producer_merge.rs` | The combination the first three could not measure — a block-filling pileup feeding a merge that skips 99 positions in 100. Also found that a four-byte summary number per locus captures 90 % of that benefit without a columnar merge. |

**The one to read first if the deferred decision comes back** is
`sketch4_columnar_producer_plus_fold.md`: it carries the working form of the change (blocks
inside the pileup, a four-byte summary column, an otherwise conventional merge) and the warning
that a block without that summary column makes the merge *worse* than today.

## `generic_walk_performance/` — the measurements behind the shipped pileup changes

Nine reports from the structural performance work of 2026-08-05, whose conclusions are
[../../../reports/reviews/perf_ng-generic-walk_2026-08-05.md](../../../reports/reviews/perf_ng-generic-walk_2026-08-05.md).
The code they describe is committed; these are the measurements, including the ones that
refuted a change.

| file | what it measured |
|---|---|
| `census.md` | What one read at one position costs, and how often the walk's general machinery is needed at all. |
| `common_column.md` | Answering an ordinary covered base in scalars instead of through the general fold. |
| `mate_overlap_sort.md` | Skipping the per-position sort when no mate pair is held. |
| `ordered_active_set.md` | Restoring the read order the walk was throwing away — and reproducing, deliberately, the 2026-05-12 regression it explains. |
| `composed.md`, `composed_full.md` | Whether those changes still add up in one tree. |
| `depth_cap.md` | Capping evidence per position instead of refusing reads at the door, after the owner called the old rule wrong. |
| `landing.md` | The assembled state and the commit sequence it was split into. |
| `price.md` | The same cost census re-taken on the shipped result — the current picture. |

Raw profiles, heap dumps and gate outputs were **not** kept; they were hundreds of megabytes and
are reproducible from the commands each report records.

---

## `str_outlier_share/` — fitting the repeat-tract outlier share (program in progress)

Scripts behind [../str_outlier_share_phase1_2026-10-18.md](../str_outlier_share_phase1_2026-10-18.md)
and [../str_outlier_share_phase2_2026-10-18.md](../str_outlier_share_phase2_2026-10-18.md), under
the plan [../str_outlier_share_plan.md](../str_outlier_share_plan.md).

| files | what it measured |
|---|---|
| `phase1.py`, `phase1b.py` | Where reads land at HG002 truth-homozygous tracts, from `pvc dev inspect-psp` dumps of psps generated with the held edge-deletion fix; one rate, a beta-binomial and two or three tract classes fitted to the per-tract counts. Paths point at the `str-group-b` worktree's `tmp/`. |
| `sim.py`, `phase2_grid_seed1.txt` | A simulator of repeat-tract reads with stray reads and noisy tracts, and a truth-free fit of the stray share; the grid over 1/10/100 samples and 3x–300x, seed 1. |
