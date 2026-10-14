# Kimura prompt: how many repeat tracts the guard leaves out of tomato25, before and after the fix

*For the owner to paste into a Claude Code session on kimura, in `/home/joxi/pop_var_caller/`. Written
2026-10-14 on branch `str-guard-pooled`.*

---

You are working in `/home/joxi/pop_var_caller/` on kimura. Read `CLAUDE.md` first. Kimura has no
container runtime: build with `cargo` directly. Scratch goes under this checkout's `tmp/` and nowhere
else. **Do not delete, move or overwrite any existing file or directory** — not psps, not earlier
results, not other worktrees. If something is in the way, stop and report it.

## What this measures

The parameters fit leaves a repeat tract out of its slippage fit when too many of the tract's reads
differ from the reference length by a number of bases that is not a whole number of repeat units.
That check is called *the guard*. The spec says it is a test of the tract: more than one such read in
ten of the reads that differ from the reference length, over the whole cohort. The code made the test
once for every sample's every read group, and left the tract out if any one failed. At about three
reads a sample, a single stray read in one sample of 2,169 was enough to remove the tract for
everyone.

On the 2026-10-13 run on this machine (branch `str-fit-c`, the tomato25 cohort, 2,169 samples) the
guard left out **45,432 of 86,688 tracts**, and only 34 of 141 repeat-length classes had enough
tracts to fit. The fix sums the counts over the cohort and tests once. We want the same counts with
the fix.

## Which branch

The tomato25 psps were written in psp format 1, which `str-fit-c` reads and `main` no longer does.
So use **`str-guard-pooled-fit-c`**: `str-fit-c` at `e03017c9` with the fix and a counting tool on
top. (`str-guard-pooled` is the same fix on `main`; it cannot read these psps.)

## Steps

1. `git fetch origin`, then make a worktree for the branch under `tmp/`, as the earlier run did with
   `tmp/pvc_str_fit_c`:

   ```
   git worktree add tmp/pvc_str_guard origin/str-guard-pooled-fit-c
   ```

   If `tmp/pvc_str_guard` already exists, stop and report; do not remove it.

2. Build the counting tool, in that worktree:

   ```
   cd tmp/pvc_str_guard && cargo build --release --example ng_str_guard_census
   ```

3. Read `tmp/str_fit_c_kimura_span9.sh` (in the main checkout) for the `--reference` and `--psp`
   arguments the 2026-10-13 counting fit used. It gave no `--catalog`; give none here either. Check
   free memory (`free -g`) and that no other large job is running. The tool reads every sample's
   repeat-tract evidence once; it fits nothing.

4. Run it with `SKIP_GATHER=1`, which stops before arranging the evidence for a fit (that step held
   19.4 GiB on this cohort and is not needed for the counts). `--output` is required by the shared
   flags and never written; give a path that does not exist:

   ```
   SKIP_GATHER=1 /usr/bin/time -v tmp/pvc_str_guard/target/release/examples/ng_str_guard_census \
       --reference <as in the script> --psp <as in the script> \
       --output tmp/str_guard_unused.toml \
       > tmp/str_guard_kimura.txt 2> tmp/str_guard_kimura.log
   ```

   Do not set `FIT_CLASSES`; it fits classes and takes hours on this cohort.

## What to report

- The whole of `tmp/str_guard_kimura.txt`. Its two count lines are the before (each sample's read
  group tested on its own) and after (tested once over the tract): tracts left out by the guard,
  tracts kept, classes with at least 8 tracts with reads, and (tract, sample) pairs with reads.
- **The before line must say 45,432 tracts left out** of 86,688. If it does not, say so first: the
  tool and the earlier run are then reading different evidence, and the after numbers mean nothing
  until that is explained.
- The three distribution lines: how the share of non-whole reads is spread among tracts kept by both
  tests, tracts the fix brings back, and tracts still left out.
- Wall time and peak memory, from the `/usr/bin/time -v` lines at the end of the log.
- The commit the worktree was built at (`git -C tmp/pvc_str_guard log --oneline -1`).

Write the report as `tmp/str_guard_kimura_report.md` and print it. Commit nothing.
