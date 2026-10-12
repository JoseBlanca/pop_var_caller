# Allele-length shapes at ±9 on tomato2, with the edges free (track C, 2026-10-12)

**Cohort:** tomato2, 59 tomato accessions over 32 Mb, about 3 reads a position. Every figure below
is from this cohort. It holds 45 repeat classes with enough tracts to fit; 12 have 100 tracts or
more, and those are all one-base runs of 8 to 13 repeats, two-base runs of 6 to 9, and the
three-base class at 6. **It cannot answer the question for long one-base runs**: those classes hold
16 to 60 tracts each.

Branch `str-fit-c`. Research plan: `doc/devel/ng/research/str_fit_research_plan.md`, avenue C.

## What was asked

A **repeat class** is every tract with one motif length and one reference repeat count. The fit
estimates the share of the class's chromosomes at each allele length, as offsets −9..+9 repeats
from the reference length. Reads are recorded only to ±8, so the two outermost offsets are
catch-alls: −9 means "9 or more repeats shorter", +9 "9 or more longer". The question is whether a
shape with a few numbers can replace those 19 free shares.

On kimura at ±6 (2026-10-11), the six-number shape `ref_step_power6` was adequate up to about 11
repeats and failed on long one-base runs, where the free share at the −6 edge was larger than at
−5: a bump no falling tail can follow.

## What changed in the shape script

`scripts/research/str_allele_length_shapes.py` now offers every shape in two variants:

- **plain**, as before: the tail runs to the edge and treats it as one more length;
- **`_edges`**: the shape covers −8..+8, and the shares at "−9 or further" and "+9 or further" are
  two more free numbers outside the tail. `ref_step_power6_edges` has eight numbers.

`--shapes plain` runs the old shapes alone. On two classes of the earlier ±6 rows, the plain shapes
reproduce the old script's log-likelihood losses and changed calls exactly.

The two measures mean what the 2026-10-11 report says. The **log-likelihood loss** is the free
fit's composite log-likelihood minus the shape's, over the same rows; the edge allele is scored at
exactly ±9, as the genotype rows give it. **Calls changed** counts rows (a sample at a tract) whose
most likely genotype differs between the free shares and the shape's; a genotype with an edge allele
is one call. The script also now measures the **free fit's own noise** on the same terms: the mean
calls changed when the free fit is swapped for one of its resampled refits (tracts resampled with
replacement). A shape that changes fewer calls than that is within the noise.

Further options, for the kimura run: `--jobs` (classes in parallel), `--max-tracts`, `--numbers`
(each class's fitted numbers as shares and powers), `--shapes none` (the free fit and its noise
only).

## Results

Shapes fitted on at most 1,500 tracts a class, with 20 resamples. Each cell is calls changed per
10,000 rows, with the median log-likelihood loss over the group's classes in brackets.

| group | classes | rows | `ref_power4` | `ref_power4_edges` | `ref_step_power6` | `ref_step_power6_edges` | free fit's noise |
|---|---|---|---|---|---|---|---|
| ≥100 tracts | 12 | 373,315 | 7.2 (70.9) | 4.6 (47.9) | 6.5 (56.0) | **3.9** (31.4) | 14.0 |
| <100 tracts | 33 | 31,555 | 180.3 (19.6) | 129.0 (10.6) | 141.3 (15.0) | **107.7** (8.3) | 272.9 |
| all | 45 | 404,870 | 20.7 (24.2) | 14.3 (15.5) | 17.0 (17.4) | **12.0** (9.9) | 34.2 |
| 1-base, 8–11 repeats | 4 | 283,624 | 6.2 | 2.9 | 6.2 | 2.9 | 2.0 |
| 1-base, 12–16 | 5 | 46,231 | 48.0 | 45.9 | 36.8 | 35.0 | 135.5 |
| 1-base, 17–21 | 5 | 4,590 | 254.9 | 268.0 | 265.8 | 259.3 | 427.9 |
| 2-base, 6–10 | 5 | 47,177 | 6.8 | 4.9 | 5.5 | 3.0 | 15.0 |
| 2-base, 11–15 | 5 | 6,350 | 40.9 | 36.2 | 67.7 | 50.4 | 287.2 |
| 2-base, 16–27 | 11 | 2,697 | 882.5 | 356.0 | 508.0 | **237.3** | 669.6 |
| 3- to 6-base | 10 | 14,201 | 18.3 | 12.7 | 9.9 | 7.7 | 48.5 |

The other shapes, over all 45 classes: `smm` 65.4, `ref_geo4` 36.6, `ref_power3` 29.5; with edges
54.4, 28.6 and 21.4.

**1. With the edges free, the shape stays within the free fit's noise in every group but one.**
`ref_step_power6_edges` changes 12.0 calls in 10,000 over all classes, against 17.0 for the plain
shape and 34.2 for the noise. The exception is one-base runs of 8 to 11 repeats, the four largest
classes: 2.9 against a noise of 2.0. That is driven by the 8-repeat class, where the edged shape
changes 53 of 81,933 calls and the plain one 14, although the edged shape's loss is lower (25.9
against 75.2). I did not find out why.

**2. The edges matter in long dinucleotides.** In two-base runs of 16 to 27 repeats, freeing the
edges halves the changed calls, from 508 to 237 in 10,000 (noise 670), and the median loss falls
from 15.0 to 6.3. In 2-base, 25 repeats the loss falls from 65.5 to 4.1. These are classes of 11 to
34 tracts.

**3. Long homopolymers: no answer from tomato2.** One-base runs of 17 to 21 repeats change about 260
calls in 10,000 with or without the edges, against a noise of 428. Five classes of 16 to 32 tracts
cannot tell the shapes apart.

**4. The shape still loses likelihood where the data are rich.** At the 1,500-tract classes the
edged shape's loss is 26 to 48 log-likelihood units, so the shape is distinguishable from the free
shares there even though the calls barely move.

### The edge shares across repeat count

`ref_step_power6_edges`, the share at "9 or more shorter" (the shape's, then the free fit's):

| class | tracts | shape | free fit |
|---|---|---|---|
| 1-base, 8–15 | 88–1,500 | 0.0000–0.0002 | 0.0000–0.0012 |
| 1-base, 16 | 60 | 0.005 | 0.007 |
| 1-base, 17 | 32 | 0.000 | 0.007 |
| 1-base, 18 | 22 | 0.016 | 0.016 |
| 1-base, 19 | 30 | 0.022 | 0.000 |
| 1-base, 20 | 23 | 0.006 | 0.006 |
| 1-base, 21 | 16 | 0.012 | 0.015 |
| 2-base, 6–15 | 20–426 | 0.0000 | 0.0000 |
| 2-base, 16 | 34 | 0.037 | 0.028 |
| 2-base, 17 | 25 | 0.015 | 0.022 |
| 2-base, 18 | 25 | 0.043 | 0.044 |
| 2-base, 19 | 19 | 0.032 | 0.036 |
| 2-base, 20 | 21 | 0.081 | 0.072 |
| 2-base, 21 | 13 | 0.215 | 0.114 |
| 2-base, 22 | 17 | 0.015 | 0.000 |
| 2-base, 23 | 16 | 0.230 | 0.127 |
| 2-base, 24 | 11 | 0.161 | 0.180 |
| 2-base, 25 | 21 | 0.384 | 0.397 |
| 2-base, 27 | 12 | 0.369 | 0.196 |

- **Two-base runs:** zero up to 15 repeats; then it rises with repeat count, from about 0.03 at 16
  to 19 to 0.16 to 0.38 at 23 to 27. In the longest dinucleotides, a large share of the cohort's
  alleles sit 9 or more repeats shorter than the reference. The trend is clear; single classes
  wander around it (0.015 at 22 repeats), as classes of 11 to 21 tracts will.
- **One-base runs:** zero up to 15 repeats; 0.005 to 0.022 from 16 to 21, with no trend visible in
  classes this thin. At 17 and 19 repeats the shape and the free fit disagree (0.000 against 0.007;
  0.022 against 0.000).
- **"9 or more longer"** stays below 0.02 in every class.

**The shorter tail's power** (how fast the shares fall beyond two repeats short) no longer drops to
zero in long one-base runs once the edge is free: 0.53, 4.6, 6.5, 3.7 and 0.39 at 17 to 21 repeats,
where the plain shape gives 0.82, 2.8, 0.00, 2.4 and 0.00. It wanders. In two-base runs of 21
repeats and more it is still zero with the edges free: the shares inside −8..−1 are flat, so those
classes spread across many shorter lengths, not only past the edge.

### A shape that ignores impossible lengths

A tract of 8 repeats cannot carry an allele 9 or more repeats shorter, and a 6-repeat dinucleotide
none 7 or 8 shorter. The shapes put a little mass there anyway (`ref_step_power6_edges` gives 0.0005
at +9 and 0.0000 at −9 for 1-base, 8 repeats; the tail itself runs to −8). The free fit puts
almost none. The effect on calls is small here, but a shape for the caller should end its shorter
tail at the tract's own length.

## Cost

Measured on the Mac (Apple silicon) while another job held the load average at 45 to 65, so the
wall times are upper bounds.

| step | wall time | peak memory | notes |
|---|---|---|---|
| counting fit, `--pairs`, ±9 | 62.7 s | 1.1 GB | 45 classes; one spread at the bound, 3-base at 13 repeats (9 tracts) |
| genotype rows, ±9 | 6.1 s | 1.3 GB | 1.1 GB on disk, against 533 MB at ±6 |
| all ten shapes, 20 resamples, 8 processes | 1 h 45 min | 1.4 GB | memory as `time -l` reports it, the largest single process |
| the noise alone (`--shapes none`), 6 processes | 12 min 15 s | 1.4 GB | same |

## Uncertain, or not done

- The edged shape changes more calls than the plain one at 1-base, 8 repeats (53 against 14),
  despite fitting better. Not investigated; Nelder–Mead in eight numbers may stop short, or the
  extra calls may sit at near-ties.
- The per-class random seed changed (it was one stream for all classes), so the bootstrap counts
  ("shares within 2 bootstrap errors") differ slightly from earlier runs on the same rows. Losses and
  changed calls do not depend on it except through the 1,500-tract draw.
- tomato2 has no well-measured long class. The kimura run (`kimura_prompt_shapes_span9.md`) is the
  test.

## Files

On this Mac, in the worktree's `tmp/` (not committed): `tomato2_span9.sh` (the commands),
`tomato2_counted_span9.tsv` and `.log`, `rows_tomato2_span9/`, `tomato2_shapes_span9.txt`,
`.numbers.tsv` and `.log`, `tomato2_noise_span9.txt`, `summarise_span9.py` (the tables above).
