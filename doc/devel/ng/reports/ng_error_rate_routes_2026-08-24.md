# ng — the two error-rate routes on one cohort, and the defect that stopped the comparison short

**2026-08-24**, branch `ng-error-rate-routes`. The comparison
[`parameter_prepass.md`](../spec/parameter_prepass.md) §4.1 has left open since the second route was
built, and [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §8 sets out as
its third measurement — both routes run on the same real reads.

---

## 1. What this establishes, and what it does not

**Pooled over the whole cohort the two routes agree to 2 parts in 1,000. Per library the comparison
could not be made at all, because the cohort tool merges all 63 libraries into one.**

| | marginal error rate | Phred |
|---|---|---|
| cohort route, all 63 libraries pooled into one | 3.9813 × 10⁻³ | 24.00 |
| per-sample route, read-weighted mean over 63 | 3.9721 × 10⁻³ | 24.01 |

The cohort tool prints its two rates to five and four decimal places, so its marginal is good to
about 1 part in 700 — **the two agree to within what this run can resolve.**

**§4.1's rule is that agreement to within rounding settles the question. This does not settle it**,
for the reason in §3 below: the number on the left is sixty-three libraries averaged into one before
the fit ever ran, and the libraries' own rates span 15.6-fold.

## 2. The defect that stopped it

**A read-group identifier is minted per alignment file, so a cohort assembled one file at a time
gives every single-library sample the same identifier.** All 63 tomato accessions arrive claiming
index 0, and `CohortCensusEvidence::read_groups`
([`joint/census.rs:1414`](../../../../src/ng/parameter_estimation/joint/census.rs)) unions them into
a set. Sixty-three libraries become one, with no duplicate check and no refusal, and the fit prints
one line where it should print sixty-three:

```
read group ReadGroupId(0): a read misreads at 0.00334 at an ordinary position and 0.0237 at a mismapped one
```

**It is in shipped code, not only in the measurement tool.**
[`estimate_contamination.rs:674`](../../../../src/pop_var_caller_exp/estimate_contamination.rs) makes
the same one-file call as
[`ng_joint_records_walk.rs:1007`](../../../../examples/ng_joint_records_walk.rs).

**Nothing wrong reaches a file today.** The fit's contamination output is keyed by sample *and*
library so it survives intact; only the error rates, keyed by library alone, collapse — and
`estimate-contamination` emits contamination and discards the rates. What it costs is that **the
cohort route's per-library error rate has never been produced by any driver that exists**, which is
the number the read likelihood would consume if that route is the one kept.

The fix is handed off separately: a cohort must refuse identifiers that collide across samples rather
than merging them, and the two callers must mint once over every alignment. A read group belongs to
exactly one sample — `build_read_groups` mints one entry per declared read group per file and
`group_by_sample` files each under the single sample its header names — so a duplicate is always the
defect and there is no legitimate case to carve out.

## 3. Why the pooled agreement is weaker evidence than it looks

**The two sides of the table are not the same operation.** The left is *one fit over sixty-three
samples' pooled evidence*. The right is *sixty-three separate fits, averaged afterwards, weighted by
reads*. A pooled fit on heterogeneous data is not in general the average of the separate fits, so
the 2-parts-in-1,000 gap mixes two things that cannot be told apart from it:

- how much the two **routes** differ, which is the question; and
- how much **pooling** differs from averaging on this cohort, which is not.

**And the cohort is heterogeneous**: the per-sample route puts the libraries' own rates 15.6-fold
apart, 6.0 in ten thousand to 9.3 in a thousand
([`ng_error_rate_spread_2026-08-24.md`](ng_error_rate_spread_2026-08-24.md)). Two routes can agree
on a cohort mean and disagree completely about which library is which, and that is exactly the
disagreement a per-library rate exists to carry.

**So the honest reading is: nothing so far suggests the two routes disagree about the cohort as a
whole, and the question §4.1 asks is untouched.**

## 4. What each run was

Both walked the same 8.0 Mb of BED over the same 63 CRAMs, against the same reference and the same
repeat catalog.

**Per-sample route** — [`examples/ng_histogram_error_rates.rs`](../../../../examples/ng_histogram_error_rates.rs),
new with this work; the route had no real-data driver taking more than one sample before it. One
walk per accession over every site in the BED, fitted by the pre-pass's own entry point
`estimate_generic_parameters`. **5,457,823,287 read-positions** in total. All 63 fitted from their
own sites, none borrowed, none railed at an end of the error-rate ladder.

**Cohort route** — [`examples/ng_joint_records_walk.rs`](../../../../examples/ng_joint_records_walk.rs),
unchanged. **1,999,404 selected positions**, the same set in every sample — one position in four of
the BED, repeat tracts excluded — against a target of 2,000,000. The fit converged in **30 passes**,
log-likelihood −43,405,858, **926.9 s**; the walk that fed it took about forty minutes and held 2.4
to 3.0 MB of census records per sample. It also reported a cohort mismapped-position share of 0.0315
and an expected heterozygosity of 4.155 per kilobase.

**Three ways the two runs are not like-for-like**, none of them fixable by running them again:

- **Site sets.** Every site in the BED against one position in four, chosen, with repeat tracts
  excluded.
- **Depth ladders.** The two bin depth differently above 8 reads a position, which bites the deeper
  accessions — this cohort runs 2.5× to 28.6×.
- **Inbreeding.** The per-sample route runs with the coefficient supplied at zero, because the runs
  model needs 3,000 separate 100 kb windows and this BED touches about 80; the cohort route fits its
  homozygote excess freely.

That they land 2 parts in 1,000 apart *despite* those three differences is the strongest thing in
this report. It is also why the per-library comparison is worth the wait rather than being
approximated from this one.

## 5. What the comparison should be judged on when it can run

**In 39 of the 63 accessions the per-sample route refused its second class of site** — the sample
asked for a noise level outside the range the model can represent, so it was fitted as though every
site were ordinary. Four more ran the coupled fit out of iterations
([`ng_error_rate_spread_2026-08-24.md`](ng_error_rate_spread_2026-08-24.md) §4).

**The cohort route holds its mismapped share for the whole cohort rather than per sample**, which is
precisely the mechanism that should rescue a sample too shallow to fit its own. Whether it does is
the sharpest question these two routes can be asked, and it is a per-library question. It cannot be
asked at all until the collision is fixed.

## 6. What is owed

- **The per-library comparison**, blocked on the read-group fix.
- **The other two axes §8 asks for** — memory and wall clock at realistic scale — not measured here.
- **The other cohort.** Everything above is one crop, 63 samples, 2.5× to 28.6×. The human benchmark
  at one sample and at 300× is the other end of both axes this caller commits to, and neither route
  has been run there for this purpose.
- **The recorded run of 2026-08-13** (`joint_records_on_real_alignments_2026-08-13.md`) reports a
  per-read-group noise line from this same tool. It has the same collision, and its line should be
  read as the cohort pooled.
