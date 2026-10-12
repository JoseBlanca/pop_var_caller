You are on kimura, in the pop_var_caller checkout. Read `CLAUDE.md` first. kimura has no container
runtime: build and run with `cargo` directly, and write every scratch file under the checkout's
`tmp/`. Python runs as `uv run --no-project --with numpy --with scipy python ...` (or plain
`python3` if numpy and scipy are installed).

## What this run is for

The owner is looking for a few-number description of how common each allele length is in a
**repeat class**: every repeat tract with one motif length and one reference repeat count. The fit
estimates, for each class, the share of the class's chromosomes at each allele length, as offsets
from the reference length. That window has been **−9..+9 repeats since 2026-10-12** (it was −6..+6).
Reads are recorded only to ±8, so **the two outermost classes are catch-alls**: −9 holds every
allele 9 or more repeats shorter, +9 every allele 9 or more longer. That makes 19 free shares a
class. The question is whether a shape with a few numbers can replace them.

**What this machine measured on 2026-10-11, at ±6**, on the tomato25 cohort (2,169 accessions,
about 3 reads a position), in your `tmp/str_shapes_kimura_2026-10-11.md`:

- A six-number shape, `ref_step_power6`, was adequate up to about 11 repeats. Its numbers: the
  share at the reference length; how the rest splits between shorter and longer; on each side, the
  weight two repeats away relative to one away, and a power for how the weights fall beyond that.
  Swapping the free shares for it changed 2 to 6 single-sample calls in 10,000 for one-base runs of
  8 to 11 repeats, no more than the free fit's own noise.
- For one-base runs of 17 to 23 repeats it changed 37 to 123 calls in 10,000, 1.3 to 1.8 times the
  noise. The free share at −6 was larger than at −5 from 15 repeats up (0.114 against 0.021 at 23
  repeats). A tail that falls with distance cannot follow a bump at the edge, so the fit flattened
  the whole shorter tail. The report put this down to the ±6 edge collecting every allele further
  away.

**This run repeats that measurement at ±9, with a new variant of each shape that gives the two edge
classes their own free numbers.** The variant's tail covers −8..+8; the shares at "−9 or further"
and "+9 or further" are two more numbers, outside the tail. So `ref_step_power6_edges` has eight
numbers.

On tomato2 (59 accessions, 32 Mb, about 3 reads a position) the edge variant changed 12.0 calls in
10,000 over all 45 classes, against 17.0 for the plain shape and 34.2 for the free fit's own noise;
in long dinucleotides (16 to 27 repeats) the edge shares came out between 0.015 and 0.38 and rose
with repeat count. See `doc/devel/ng/research/str_fit/trackC_tomato2_span9_2026-10-12.md` in this
branch. tomato2's homopolymer classes above 15 repeats hold 16 to 60 tracts each, too few to decide
anything there; kimura has many more.

## The questions

(a) **Is the edge-free shape adequate, class by class, at ±9?** For `ref_step_power6_edges`, and
    `ref_step_power6` beside it: one-sample calls changed per 10,000 rows, and the composite
    log-likelihood loss, by class; grouped by tract count (100 or more tracts, fewer than 100), as
    in the 2026-10-11 table. Set the changed calls against the free fit's own noise, as the
    2026-10-11 follow-up did. The script now measures that noise itself, on the same tracts and
    resamples (the line `free fit vs its own refits (noise)` and the per-class `noise` column), so
    `tmp/str_noise_baseline.py` is not needed.

(b) **Do the shape's numbers, including the two edge shares, follow smooth curves across repeat
    count, within each motif length?** That is what would let thin classes borrow from their
    neighbours. For each number, say whether it rises, falls or wanders, and from what to what.
    In particular: once the edge has its own share, does the power below stop collapsing to zero
    from 19 repeats up?

(c) **The counting fit's cost at ±9 on 2,169 samples.** Wall time and peak memory. At ±6 it took
    15 min 38 s and 31.5 GiB peak (33,034,756 kB maximum resident). Say whether either grew, and
    by how much.

## The run

### 1. The code

```
git fetch origin
git worktree add ../pvc_str_fit_c origin/str-fit-c     # or, if worktrees are not allowed here:
                                                         # git worktree add tmp/pvc_str_fit_c origin/str-fit-c
cd <that worktree>
git log -1 --format='%h %s'                             # note the commit in the report
cargo build --release --example ng_str_counting_fit --example ng_str_genotype_rows
```

On 2026-10-11 the project rules ruled out `../`, and the worktree went to `tmp/pvc_str_shapes`; do
the same if that still holds. Run every command below from inside the worktree, and point every
output at the **main checkout's** `tmp/` (written `$T` below), so the results sit beside the
2026-10-11 ones.

### 2. The inputs (fill in from the 2026-10-11 run)

Take them from `tmp/str_shapes_kimura.sh` in the main checkout, which that run used:

- `REF`: its `--reference`
- `CATALOG`: its `--catalog`
- `PSPS`: its `--psp`
- `PARAMS`: its `--fitted`, which was `tmp/fit_precision_kimura_2026-10-06.parameters.toml` (each
  sample's homozygote excess)

Write them into a launcher, `$T/str_fit_c_kimura_span9.sh`, so the report can name it.

### 3. The counting fit, every class, at ±9

```
/usr/bin/time -v ./target/release/examples/ng_str_counting_fit \
    --reference REF --catalog CATALOG --psp PSPS --output $T/unused --fitted PARAMS \
    --pairs --allele-span 9 \
    > $T/str_counted_kimura_span9.tsv 2> $T/str_counted_kimura_span9.log
```

`--allele-span 9` is the default now; it is written out so the log says it. If the program refuses
the psps ("another selection", or a psp that does not decode), stop and report that, with the
message. Do not regenerate the psps.

From `time -v`: "Elapsed (wall clock) time" and "Maximum resident set size". Also list any class
whose `concentration` sits at a bound (`0.001` or `10000`; on tomato2 one class did, 3-base at 13
repeats).

### 4. The genotype rows at ±9

```
./target/release/examples/ng_str_genotype_rows \
    --reference REF --catalog CATALOG --psp PSPS --output $T/unused --fitted PARAMS \
    --allele-span 9 --classes $T/str_counted_kimura_span9.tsv --out-dir $T/str_rows_kimura_span9 \
    2> $T/str_rows_kimura_span9.log
```

**Disk:** at ±6 these rows took 46 GB. A row at ±9 holds 190 genotypes against 91 at ±6, so expect
about 95 GB. Check `df -h` on the checkout first; if there is not room, delete
`tmp/str_rows_kimura` (the ±6 rows) only after asking the owner, and otherwise stop and report.
Note the final size (`du -sh`).

### 5. The shapes

The script is `scripts/research/str_allele_length_shapes.py` in this branch. It now runs classes in
parallel itself (`--jobs`), caps tracts per class itself (`--max-tracts`), and writes each class's
fitted numbers (`--numbers`), so the 2026-10-11 helper copies under `tmp/` are not needed.

**Measure the cost on one class before the full run.** The largest class (one-base, 8 repeats) took
3 h 34 min at ±6 with 300 tracts, 10 resamples and five shapes. Rows at ±9 are about twice as wide,
and each edge variant has two more numbers to fit. Make a directory of links holding
`genotypes.tsv` and the one file `1_8.bin`, and run:

```
mkdir $T/str_rows_one_class
ln -s $T/str_rows_kimura_span9/genotypes.tsv $T/str_rows_kimura_span9/1_8.bin $T/str_rows_one_class/
/usr/bin/time -v uv run --no-project --with numpy --with scipy python \
    scripts/research/str_allele_length_shapes.py $T/str_rows_one_class 0 10 \
    --max-tracts 300 --shapes ref_step_power6,ref_step_power6_edges,ref_power4_edges \
    > $T/str_shapes_kimura_span9_1_8.txt 2> $T/str_shapes_kimura_span9_1_8.log
```

If that class takes more than about 6 hours, report its time and ask the owner before starting the
full run. Otherwise run all classes, 16 at a time (set `OMP_NUM_THREADS=1` so numpy does not start
threads of its own in each process):

```
OMP_NUM_THREADS=1 /usr/bin/time -v uv run --no-project --with numpy --with scipy python \
    scripts/research/str_allele_length_shapes.py $T/str_rows_kimura_span9 0 10 \
    --max-tracts 300 --jobs 16 \
    --shapes ref_step_power6,ref_step_power6_edges,ref_power4_edges \
    --numbers $T/str_shapes_kimura_span9.numbers.tsv \
    > $T/str_shapes_kimura_span9.txt 2> $T/str_shapes_kimura_span9.log
```

The script prints each class to the log as it finishes. Its output has three parts: the free shares
by class (the two ends are "this far or further"); one line per shape over every class, then the
free fit's noise on the same terms; and per class, the noise (mean calls changed when the free fit
is swapped for one of its resampled refits) and, for each shape,
`log-likelihood loss / largest share miss / rows whose call changes`. The
numbers file has, per class and shape, the fitted numbers as shares and powers
(`ref_share`, `shorter_split`, `two_step_weight_below`, `power_below`, `two_step_weight_above`,
`power_above`), the shape's two edge shares (`edge_below`, `edge_above`), and the free fit's edge
shares beside them (`free_edge_below`, `free_edge_above`). **In an `_edges` shape, `ref_share` and
`shorter_split` are shares of what the edges leave**, not of the whole; `inner_ref_share` is the
same for plain and edged shapes and compares them. A power printed as about 5.2e21 sits at the
script's cap: the side it describes holds almost nothing, so the power is not measured.

What the two measures mean (unchanged from 2026-10-11):

- **Composite log-likelihood loss:** the free fit's composite log-likelihood minus the shape's, over
  the same rows. It depends heavily on which 300 tracts are drawn, so read its trend, not single
  values.
- **One-sample calls changed:** rows (a sample at a tract) whose most likely genotype under the free
  shares differs from the one under the shape. A genotype with an edge allele is one call, whatever
  distance past the edge the allele really sits at.

If you need the old shapes too (`smm`, `ref_geo4`, `ref_power3`, `ref_power4`), `--shapes all` runs
all ten; it costs about three times as long as the three above.

## The report

Write `$T/str_shapes_kimura_span9_2026-10-12.md` (use the actual date) and print its path. The
owner brings it back. Lead with the answers to (a), (b) and (c), with the numbers that decide each.

1. **(a)** The table of 2026-10-11, now with `ref_step_power6` and `ref_step_power6_edges` (and
   `ref_power4_edges`) side by side: changed calls per 10,000, over classes with 100 or more tracts
   and over the rest, with rows and classes counted in each group; the worst class of each. Then
   the change rate by class for one-base runs, as the 2026-10-11 list gave it (8–11, 12–16, 17–23
   repeats), and for two-base runs. Then the noise comparison.
2. **(b)** For `ref_step_power6_edges`, a table of its numbers by motif length and repeat count, for
   classes with 100 or more tracts, with the free fit's edge shares beside the shape's. Say how each
   number moves with repeat count. Compare the power below with the 2026-10-11 table, where it was
   0.000 from 19 repeats up.
3. **(c)** The counting fit's wall time and peak memory at ±9, against 15 min 38 s and 31.5 GiB at ±6;
   the classes fitted; any concentration at a bound. The genotype rows' size and time. The shapes
   run's time and the largest class's time.
4. Every path: launcher, logs, tables, the numbers file, the commit you ran.

Do not change committed code and do not commit anything. Anything you write goes in `tmp/`.
