# C4 — what calling repeat tracts is worth, and the one thing it broke

**Date:** 2026-09-02. **Plan:** [`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md)
Milestone C step C4. **Runs:** `benchmarks/giab/src/run_ng_per_sample.sh` at 30× and 50×, all
three GIAB per-sample callsets, `--defaults`, scored with
`benchmarks/giab/src/score_ng_recall.sh`.

---

## The headline

**Indel recall goes from 0.673 to 0.915 at 30× and from 0.676 to 0.939 at 50×**, pooled over the
three samples. SNP recall goes from 0.974 to 0.980 and from 0.979 to 0.984. Nothing about the
reads or the model changed — what changed is that a repeat tract is now called.

**And the same variant is written twice.** Indel precision falls from 0.987 to 0.816 at 30×, and
**62 of the 68 new false calls are a record the file already holds**: the SNP/indel path emits an
insertion anchored at the base beside a tract, and the tract path emits the same insertion as a
length change. Counted once, precision is 0.981 against the baseline's 0.987. The defect is the
duplication, not the tract calls.

## The numbers

Pooled over HG002, HG003 and HG004, each on its own 100-interval confident BED against its own
GIAB truth VCF. **Both columns are from the same script and the same inputs**; the baseline was
re-run today from the pre-dispatch binary so that the two depths are the same vintage.

| depth | class | | before | after | after, counting a duplicate once |
|---|---|---|---:|---:|---:|
| 30× | indels | recall | 0.673 (222/330) | **0.915** (302/330) | 0.915 |
| | | precision | 0.987 (222/225) | 0.816 (302/370) | **0.981** (302/308) |
| | SNPs | recall | 0.974 (2 008/2 061) | **0.980** (2 019/2 061) | |
| | | precision | 0.993 | 0.993 | |
| 50× | indels | recall | 0.676 (223/330) | **0.939** (310/330) | 0.939 |
| | | precision | 0.987 (223/226) | 0.809 (310/383) | **0.978** (310/317) |
| | SNPs | recall | 0.979 (2 018/2 061) | **0.984** (2 029/2 061) | |
| | | precision | 0.992 | 0.992 | |

**Against the stated bar**, which is the production caller on the same samples and regions
([the routing report](../ng_str_routing_recovery_2026-09-02.md)): SNP 0.987, indel 0.930 at 30×.
ng is now 0.980 and 0.915.

**And against ng's own best, which is where the honest gap is.** Before the region clip landed —
a record now stops at the end of the region it was walked for, which cost 90 of 312 true indels —
ng's indel recall at 30× was 0.946. Tract calling brings it back from 0.673 to 0.915, so **about
a third of what the clip cost is still missing**, and the tract path is what was expected to
return all of it.

## The duplication, measured

After `bcftools norm -m -any` — left-aligned and split, which is how the scorer compares — the
query file holds the *same* `(CHROM, POS, REF, ALT)` twice:

| depth | HG002 | HG003 | HG004 | total | new false indel calls |
|---|---:|---:|---:|---:|---:|
| 30× | 28 | 15 | 19 | **62** | 68 |
| 50× | 29 | 17 | 20 | **66** | 73 |

**Nine in ten of the new false calls are a duplicate**, and every distinct false position also
carried a true call, at a position where the truth set does have an indel. Nothing is being
invented on quiet ground.

Here is one, from HG002 at 30× — two ng records, one base apart, describing one inserted `A`:

```
chr1 4927622 T             TA              PASS  1/1   ← the SNP/indel path
chr1 4927623 AAAAAAAAAAAA  AAAAAAAAAAAAA   PASS  1/1   ← the tract path, STR/RU=A/PERIOD=1
```

Left-aligned they are the same variant. The truth set carries it once, so one of the two is a
true positive and the other is a false one.

**Where it comes from.** The routing puts the tract on the repeat path and the base immediately
before it on the generic path; the generic mint opens a record at that base for an insertion
whose content is the tract's, and the tract generator describes the same event as a length
change. Before this milestone the second record did not exist, so the duplication could not
happen.

## What is not affected

**SNP precision does not move**: 0.993 at 30× before and after, 0.992 at 50×. The duplication is
an indel phenomenon, which is what the mechanism above predicts — a tract's alleles differ from
the reference in length.

**Thread-invariance holds with tract records in the file.** The lib suite's
`the_tract_path_is_byte_identical_at_every_thread_count` sweeps 1, 2, 4, 8 and 16 threads three
times each at two merge widths, on the fixture that varies inside its tract, and its non-vacuity
anchor is now a record over the tract's own ground rather than a set-aside count that is zero.

**The dashboards need no edit.** `accuracy_dashboard.py` reads the results directories from disk,
and the new calls are in `results/per_sample/{30x,50x}/ng`; the pre-dispatch calls are kept beside
them in `ng_before_tract_calling`, from the same script and the same inputs.

## The tract outcomes the runs report

| depth | sample | built | called | notPeriodic |
|---|---|---:|---:|---:|
| 30× | HG002 | 75 | 74 | 1 |
| | HG003 | 57 | 57 | 0 |
| | HG004 | 74 | 72 | 2 |
| 50× | HG002 | 98 | 98 | 0 |
| | HG003 | 79 | 79 | 0 |
| | HG004 | 98 | 96 | 2 |

No tract was refused for `tooManyAlleles`, for a candidate carrying no whole motif copy, or as a
bundle, on this benchmark's ground.

## What this leaves for a decision

**The duplication is a design question that spans two plans and is not this step's to settle.**
Three places could fix it, and they are not equivalent:

- **the routing**, by giving a tract's region the anchor base beside it, so the generic mint never
  opens a record there — cleanest, and it belongs with typed regions;
- **the generic mint**, by refusing to open a record whose footprint reaches into a repeat region
  — narrower, and it puts repeat knowledge in the SNP/indel path;
- **the output**, by dropping a record another record already describes — cheapest and the least
  honest, since it hides which path was right.

**Recommended: the routing.** It is where the boundary is already decided, it is the only one of
the three that prevents the second record rather than removing it, and the region clip that
created the seam landed there too.

## Validation

The runs are the shipped release binary built from this branch, in the container. Scoring is
`score_ng_recall.sh`, which follows `accuracy_dashboard.py`'s method: both truth and query
restricted to the sample's BED, left-aligned and split with `bcftools norm -m -any`, filtered to
one class, intersected on position and alleles.
