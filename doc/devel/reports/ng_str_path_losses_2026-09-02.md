# Where ng loses the variants inside repeats

**Date:** 2026-09-02. Measured on the GIAB per-sample benchmark at 30× and 50×, three
samples (HG002, HG003, HG004) over their own confident regions — 2,061 truth SNPs and 330
truth indels pooled.

Follows [`ng_first_calling_benchmark_2026-09-02.md`](ng_first_calling_benchmark_2026-09-02.md),
which measured that nine in ten of the truth variants ng misses are in stretches of the
reference the caller builds no locus over. This says where in the pipeline that happens, and
how much of it is a repeat problem at all.

**Two words this report leans on.** A **locus generator** is the code that turns one stretch
of the reference plus the reads lying over it into a *locus* — a place with a list of the
distinct sequences the reads showed and how many reads showed each. It runs before any
genotyping: it builds the thing the calling loop then calls. ng has one slot for each kind of
stretch the reference is cut into, and fills only one of them. **Analysed sequence** is the
bases of the reference a run was asked to look at — the BED's spans, in bases. The code calls
an ordinary stretch — one where nothing more specific than "sequence" can be said — a
**generic** region, and that word appears in the tool output quoted below.

---

## 1. The two losses, and they are not the same size

**ng emits nothing at a repeat tract because the tract slot in its locus generator set is
empty — and it classifies about seven times more of the reference as repeat than its own
calling policy would, so seven times more sequence reaches that empty slot.**

*Sequence, measured in bases.* On HG002's confident regions the run was handed 572,037 bases
to analyse; it classified 32,577 of them as repeat and built no locus over any of them. Under
ng's own calling floors that figure is 4,930. Across the three samples the ratio is 6.6, 7.3
and 7.0. In truth variants rather than bases the ratio is smaller — 103 lost SNPs against 22,
55 lost indels against 13, about 4.5 either way — because the extra sequence the permissive
floors capture is short and variant-poor.

Those are separate faults with separate fixes:

| | what it is | what it costs, three samples pooled |
|---|---|---|
| **Routing** | the run classifies the reference with the floors the *catalog file* was stored at, not the floors ng's calling policy names | 81 of the 103 lost SNPs and 42 of the 55 lost indels sit in sequence the calling floors would leave on the generic path — which is built and working |
| **The tract path** | both repeat slots in the generator set are `Unfilled(NotImplemented)`, so a tract yields no locus | the remaining 22 SNPs and 13 indels, and every genuine repeat-length variant |

Correcting the routing alone needs no new calling code. It is a one-line change of which
criteria the run asks the catalog with, and the sequence it moves goes to a locus generator
that is built, wired in, and whose loci ng calls at 98.4% recall.

---

## 2. What the measurement is

`examples/ng_typed_region_dump.rs` prints the typed regions a run routes on — the same
catalog, the same walk, the same criteria — so each truth variant can be placed in the region
that covers it. Run twice per sample, once with each set of copy-number floors:

- **catalog floors**, `[5, 5, 4, 4, 4, 3]` over periods 1–6 with a 500 bp satellite cap. These
  are what `call-from-alignments` uses today
  ([`call_from_alignments.rs:845`](../../../src/pop_var_caller_exp/call_from_alignments.rs#L845)).
- **calling floors**, `[8, 6, 6, 6, 5, 4]` with a 100 bp cap — `SsrSegmentCriteria::default`,
  the copy counts at which a repeat starts to stutter, measured over the tomato archive on
  2026-08-10.

The probe reproduces the run exactly: on HG002's regions it puts 539,460 bases on the generic
path and 32,577 on the repeat path, which are the two figures that run's own report prints.

**The catalog's floors are storage floors, and the spec says so.**
[`repeat_catalog.md`](../ng/spec/repeat_catalog.md) §on the copy floor: the file is built below
every calling floor precisely so that "a caller can therefore move its routing floor anywhere
inside that gap by filtering, which is the question the file exists to keep open." The run does
no filtering — it routes on the file's own floors.

---

## 3. Where the truth variants sit

Three samples pooled, by the typed region covering each truth variant.

| class | floors | ordinary sequence | one repeat tract | a bundle of tracts | satellite | in a stretch with no generator |
|---|---|---:|---:|---:|---:|---:|
| SNPs | catalog (the run's) | 1 958 | 53 | 50 | 0 | **103 (5.0%)** |
| SNPs | calling | 2 039 | 15 | 4 | 3 | **22 (1.1%)** |
| indels | catalog (the run's) | 275 | 10 | 45 | 0 | **55 (16.7%)** |
| indels | calling | 317 | 6 | 5 | 2 | **13 (3.9%)** |

Two things to read off it.

**Most of the loss is in bundles, not in tracts.** A bundle is a cluster of repeats none of
which has clean flanks — 45 of the 55 lost indels and 50 of the 103 lost SNPs are in one. At
the calling floors the same sequence holds 5 and 4. Admitting more short repeats does not just
add tracts; it makes neighbours collide, and a collision swallows the sequence between them.

**Two thirds of the loss is not a repeat-length variant at all.** 103 of the 158 lost truth
variants are SNPs — substitutions that happen to sit inside or beside a repeat. They need no
repeat caller; they need to be on the generic path, where ng calls SNPs at 0.984 (§5).

## 4. The same six sites, under each set of floors

HG003, chr22. Five of the six move from a bundle nothing is built over to a generic stretch
the working caller handles, with no change to any caller:

| site | catalog floors (the run) | calling floors |
|---|---|---|
| SNP chr22:17,321,924 | `ssr_locus` 17,321,920–935 | `ssr_locus` 17,321,920–935 |
| SNP chr22:41,679,577 | `ssr_bundle` 41,679,558–595 | `generic` 41,679,257–678 |
| SNP chr22:41,680,377 | `ssr_bundle` 41,680,369–384 | `generic` 41,680,086–41,681,990 |
| SNP chr22:41,686,493 | `ssr_bundle` 41,686,472–498 | `generic` 41,686,321–632 |
| indel chr22:41,679,678 | `ssr_bundle` 41,679,670–697 | `generic` 41,679,257–678 |
| indel chr22:41,681,990 | `ssr_bundle` 41,681,983–42,007 | `generic` 41,680,086–41,681,990 |

## 5. What each caller recovers in each kind of stretch

Same three samples, same split, all three callers:

| depth | class | caller | ordinary sequence | classified as repeat |
|---|---|---|---:|---:|
| 30× | SNPs | ng | 1 926/1 958 = **0.984** | 0/103 = **0.000** |
| | | production | 1 932/1 958 = 0.987 | 102/103 = 0.990 |
| | | freebayes | 1 822/1 958 = 0.931 | 89/103 = 0.864 |
| 30× | indels | ng | 270/275 = **0.982** | 0/55 = **0.000** |
| | | production | 260/275 = 0.945 | 47/55 = 0.855 |
| | | freebayes | 256/275 = 0.931 | 45/55 = 0.818 |
| 50× | SNPs | ng | 1 937/1 958 = **0.989** | 0/103 = **0.000** |
| | | production | 1 934/1 958 = 0.988 | 102/103 = 0.990 |
| | | freebayes | 1 836/1 958 = 0.938 | 90/103 = 0.874 |
| 50× | indels | ng | 270/275 = **0.982** | 0/55 = **0.000** |
| | | production | 262/275 = 0.953 | 50/55 = 0.909 |
| | | freebayes | 258/275 = 0.938 | 47/55 = 0.855 |

**In ordinary sequence ng is the best of the three at indels and level with the production
caller at SNPs.** At 30× its indel recall there is 0.982 against the production caller's 0.945
and freebayes' 0.931; at 50× it is 0.982 against 0.953 and 0.938. That corrects the reading the
headline benchmark invites: ng's indel recall is nine to thirteen points below both *overall*,
and above both *where it calls*.

**In sequence classified as repeat it is zero, exactly, at both depths.** The other two callers are at
0.86–0.99 there. More depth changes neither figure, which is what says the gap is not evidence.

**What re-routing alone would be worth.** If the 81 SNPs and 42 indels the calling floors put
on the generic path were called at ng's own rate there, its overall recall at 30× would go from
0.935 to about 0.97 on SNPs and from 0.818 to about 0.94 on indels — against the production
caller's 0.987 and 0.930. **That is an upper bound**: those positions abut repeats and are
harder than the generic average, so the real figure is below it. It also moves 3 SNPs and 2
indels the other way, into the satellite class, which is a permanent refusal rather than an
unbuilt path — the calling floors' 100 bp cap is stricter than the catalog's 500 bp.

---

## 6. The evidence at a tract already exists

The repeat locus generator is built and it works. Run over HG003's chr22 at 30×
(`examples/ng_ssr_loci_dump.rs`, which asks the catalog with the *calling* floors), it walks
26,680 tracts and reads off tract lengths wherever
there are reads: inside the confident regions each covered tract has **18 to 29 reads that pin
its length exactly**, plus partial reads recorded separately. The three tracts with 2, 4 and 14
reads straddle a region's edge, where coverage genuinely stops.

So the loss is not evidence. At 30× a tract has ample reads, and ng already reads their
lengths; nothing takes them further.

---

## 7. The chain, stage by stage

What each stage does with a repeat tract today:

| stage | state |
|---|---|
| Region typing from the catalog | **Built.** Routes tracts, bundles and satellites correctly — with the wrong floors (§2). |
| Repeat locus generator (`SsrGenerator`, `locus_generation/ssr.rs`) | **Built and tested**, ~3,700 lines, and demonstrated above. |
| Wiring it into a run's generator set | **Missing.** Both repeat slots are `Unfilled(NotImplemented)` ([`walker.rs:1615,1617`](../../../src/ng/run/walker.rs#L1615)); only tests and examples fill them. **This is where the loci are lost** — nothing downstream ever sees a tract. **Filling this slot alone is not enough**: the calling path is unconditionally the SNP/indel one (`call_one_generic_locus`), so a tract's loci would be genotyped with the SNP/indel error model, which has no slippage term, and written as plain indel records with none of the repeat annotation the VCF encoder already knows how to produce. Wrong answers rather than a crash. |
| The merge carrying the locus kind | **Missing.** `CohortObservation` has `region`, `alleles`, `per_sample` and no kind ([`cohort_merge/build.rs:974`](../../../src/ng/run/cohort_merge/build.rs#L974)), so even with a generator there is no motif downstream, so no period, no repeat count, no ladder. Milestone A of [`candidate_alleles_ssr.md`](../ng/impl_plan/candidate_alleles_ssr.md). |
| Candidate selection at a tract | **Missing.** `src/ng/calling/allele_candidates/ssr.rs` does not exist; neither `select_ssr` nor `RepeatLadder` is written. |
| Evidence shaping | **Built** — `shape_ssr_locus`, with a note that it takes its candidate repeat counts as supplied "because the repeat-tract half of candidate selection is unwritten". |
| Read likelihood and genotyping | **Built** — the scoring context, the stutter model, the length support and the outlier weight, and the repeat arms of `summarise_condition`. |
| VCF emission | **Built** — `STR`, `RU`, `PERIOD` in INFO and `REPCN` in FORMAT, with the repeat filters already declared in the header. |
| Bundles | **No path at all.** They carry most of the loss (§3) and nothing beyond "carry them rather than delete them" has been designed. |

---

## 8. What follows

1. **Give the run a calling policy to route with**, rather than the catalog's storage floors.
   It is a change to which criteria `segments_over` asks the catalog with, and it moves about
   four in five of the lost variants onto a caller that already works. Ask the catalog with
   `StrRepeatCriteria::from(&TypedRegionConfig::default())`, which is what every development
   tool in the tree already does.
2. **Then the tract path**, in the order its own plan sets: the locus kind through the merge,
   selection at a tract, and the generator wired into the run's set.
3. **Bundles need a decision before either helps them.** Under the calling floors they are 4
   SNPs and 5 indels rather than 50 and 45, so re-routing shrinks the problem by an order of
   magnitude — but it does not answer what a caller should do with a cluster of repeats that
   has no clean flanks.

## 9. What was added

- `examples/ng_typed_region_dump.rs` — prints the typed regions a run routes on, under either
  set of copy floors, so any position can be placed in the region that covers it.
