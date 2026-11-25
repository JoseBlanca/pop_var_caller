# Refinement replay on tomato25 after the memory, solve and pipeline changes (2026-11-25)

Main f11e6897: the memory cuts, the solve speed-ups, the pipeline changes and both review rounds'
fixes. Measured on kimura with `examples/refine_replay.rs` over the kept spill of the 2,176-sample
tomato25 run, chromosome 1 at 40–47 Mb (825,782 entries), 48 threads, margins thinned at 5 records
a kb, core floor at its default (the margin floor, 500 kb). Each run is a systemd unit capped at
118 GB; "peak" is the replay's own per-window figure, "anonymous memory" the unit's sampled every
10 s. Scripts and logs: `/media/tomato25_vcfs/pvc_work/refine_scale/` on kimura (`queue_f11e.sh`,
`f11e_80/`, `f11e_65/`).

## Results

| Run | Budget | Windows (core length) | Highest window peak | Anonymous memory, max | Wall | Fill check wrong |
|---|---|---|---|---|---|---|
| 2026-11-23, unpacked, before these changes | 80 GB | 7 (0.50–0.85 Mb) | 97.6 GB | 97.5 GB | 7,587 s | 0.0800% (713 of 890,982) |
| 2026-11-23, packed, before these changes | 65 GB | stopped at 42.2 Mb | 70.9 GB | — | — | — |
| f11e6897 | 80 GB | 4 (0.85–1.45 Mb) | 103.0 GB | 102.5 GB | 5,248 s | 0.0712% (423 of 594,265) |
| f11e6897 | 65 GB | 6 (0.60–1.70 Mb) | 88.7 GB | 88.0 GB | 5,087 s | 0.0794% (551 of 694,041) |

- **The stretch that stopped the run now goes through at 65 GB.** On 2026-11-23 the smallest
  window at 42.2 Mb counted 65.9 GB; now the windows there count under the budget and run.
- **Time fell by about a third** at the same 80-GB budget (7,587 → 5,248 s), with fewer, longer
  windows.
- **Real memory per GB of budget rose.** The highest window peak was 1.29 times the budget at 80 GB
  and 1.36 times at 65 GB, against 1.16–1.22 before. The count fell by more than what is held, as
  the memory build's report predicted: at the same budget, windows grow until they hold more.
- The two budgets fit a line, peak ≈ 0.95 × budget + 27 GB. The 27 GB is held whatever the window
  size — at 48 threads, the walks' per-thread arrays, the spill reader and the allocator's slack
  are candidates; not yet broken down. Two points only: treat it as a guide between 65 and 80 GB.

## Budget for the tomato25 rerun

`snp_calling/pipeline.toml` on kimura sets 80 GB, chosen on 2026-11-23 for a peak of about 93 GB,
under 100 GB. At f11e6897 80 GB peaks at 103 GB. A budget of **70 GB** should peak at about
93 GB by the line above; 65 GB peaked at 88.7 GB. Not changed: the owner's to set.
