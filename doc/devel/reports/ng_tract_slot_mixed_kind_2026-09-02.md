# A cohort locus that mixes ordinary sequence with a repeat tract — measured, and open

*2026-09-02. Found while validating step C2 of
[`run_ssr_observations.md`](../ng/impl_plan/run_ssr_observations.md) on real data. Branch
`ng-ssr-observations`. **Milestone C is blocked on a ruling.***

---

## 1. What happens

With the repeat-tract generator wired into a run's generator set, calling HG002's 100 confident
regions at 30× **aborts**:

```
thread 'main' panicked at src/ng/run/cohort_merge/close.rs:647:
a cohort locus at ContigId(0):4916058 mixes locus kinds — Generic with
Ssr(SsrDetail { motif: Motif("GT"), … })
```

That assertion is deliberate and release-level. Its own comment gives the argument for why it
can never fire:

> It cannot happen, and structurally rather than by luck: segments are the reference's own
> partition, no observation crosses a segment boundary (run spec §4.3), so no chain of
> overlapping observations can either.

**The second clause is false**, and has been since before this branch. The SNP/indel generator
clamps an observation's **start** to its segment (`clamp.contains(locus.region.start)`,
`pileup/generator.rs`) and says nothing about where it ends. An indel's footprint spans the
anchor plus the run it covers, so a read with a deletion anchored on the last base of a generic
stretch produces an observation that reaches into whatever comes next.

Nothing noticed, because until C2 the tract slot was unfilled: there was no tract observation
for the straddling generic one to chain with.

At the position the run aborted on, the routing is:

```
chr1  4915934  4916058  generic
chr1  4916059  4916096  ssr_locus  GT  period 2  19.0 copies  purity 1.00
```

— the locus opens on the generic stretch's **last** base and reaches one further.

## 2. How often, measured

The assertion was temporarily replaced with a counted chain break and the three GIAB samples
re-run at 30×. The patch was reverted afterwards and its absence checked by grep; the benchmark
outputs it produced are set aside under `results/per_sample/30x/ng_scratch_measurement_do_not_trust/`
and are not the ones B4 scored.

| sample | loci that would mix kinds | repeat-tract loci built |
|---|---:|---:|
| HG002 | 67 | 75 |
| HG003 | 50 | 57 |
| HG004 | 60 | 74 |
| **all three** | **177** | **206** |

**So it is most tract loci, not a corner.** Roughly nine in ten of the tract loci a run builds
have a generic observation reaching into them. That is not surprising once stated: the bases
immediately beside a repeat tract are where reads carry indels.

## 3. Three ways out, and what each costs

**(a) Close the locus at the kind change** — the closer breaks the chain instead of asserting.
Two lines, and it is what the measurement above ran. **But the generic locus still *spans* into
the tract**, so the organiser's overlap rule then drops the tract locus as ground an earlier
locus already owns — the tract would be lost silently. I have not confirmed that on a fixture,
and it is the reason I am not recommending it.

**(b) Clip a generic observation's footprint at its segment's end.** This makes the assertion's
argument true instead of removing it. Cost: a deletion straddling a tract edge loses its tail,
so the SNP/indel path's calls change at every tract boundary — and B3's byte-identity pin would
fail, which is that pin doing its job on a different change rather than an objection.

**(c) The tract's ground is the tract's** — drop a generic observation that reaches into a
tract outright. Loses the whole deletion rather than its tail.

## 4. Recommendation

**(b), and measure it.** The argument is the one the merge already makes for treating the two
kinds differently: *a generic observation's span is the mapper's CIGAR taken on trust, an STR
observation's is a tract the catalog fixed before any read was looked at, and its reads were
re-aligned rather than believed*. A footprint reaching into a tract is exactly a CIGAR claim
about ground the tract path exists to re-align. Clipping keeps the part of the observation that
lies on generic ground and hands the tract's bases to the path that will look at them properly.

What it costs in calls is unmeasured and is one benchmark run away — the same 40-second GIAB
sweep B4 used, scored against the calls now on disk.

**What I need from you is which of the three**, because (b) changes SNP/indel output near every
tract and that is a design decision, not an implementation one.

## 5. Where the branch stands

C1, C2 and C3 are committed; C4's test is committed and **its checkbox deliberately not
flipped**. The unit suite is green — 5,998 passing — because no fixture in it puts an indel on
the last base before a tract. **A real run aborts**, so Milestone C is not complete and must not
be merged as it stands.

The plan's own C4 asks for byte-identical output at 1–16 threads with the slot filled, and that
now holds on a fixture where the tract path fires. It says nothing about this, because a unit
fixture had no reason to build the shape that breaks it — which is the lesson worth keeping:
**the defect was found by running the caller on a genome, not by the suite.**
