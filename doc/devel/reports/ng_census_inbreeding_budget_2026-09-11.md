# Fitting the inbreeding coefficient from the census, and what the position budget has to be

**2026-09-11**, branch `census-audit`. Programs: `examples/ng_census_inbreeding_budget.rs` (the
accuracy sweep) and `generate-psps` on the tomato benchmark (the cost), over
`src/ng/parameter_estimation/joint/census_runs.rs`. **All three were deleted the same day — see the
banner below.**

> ## ⛔ What this measured was then removed, and this report is why
>
> The owner's ruling on the strength of §3 and §5: **the census budget stays at two million and the
> runs-of-homozygosity estimator goes**, with the whole-genome histogram route it lived in
> ([`impl_plan/remove_histogram_route.md`](../ng/impl_plan/remove_histogram_route.md)). The
> reasoning is §5's last row — tripling the census moved the parameters it was already fitting by
> under 1%, so the extra 4.86 MB a sample bought only this one estimator, and the coefficient a
> caller reads is the homozygote excess instead.
>
> **So `census_runs.rs` and `examples/ng_census_inbreeding_budget.rs`, which produced everything
> below, no longer exist.** The numbers stand as the measurement the decision was taken on; the code
> is in this branch's history.

---

## 1. What this establishes

**The inbreeding coefficient can be fitted from a census, and six million positions is enough for
a typical accession of the tomato panel but not for its least heterozygous ones.** Raising the
census's generic budget from two million positions to six cuts the worst error in the coefficient
from 0.105 to 0.086 on a selfing genome at the cohort's median heterozygosity; sixteen million would
cut it to 0.026.

**What it costs, measured on four tomato accessions over the 8 Mb benchmark: 4.86 MB a sample on
disk, and — the larger number — 2.78× the peak memory and 2.13× the wall clock at the fit**, from
185.8 MB and 944 s to 516.9 MB and 2,010 s. How the fit's memory scales with cohort size is not
measured.

**And tripling the census moved the parameters it was already fitting by under 1%**, so the raise
buys nothing for those and is worth its cost only if the coefficient is going to be fitted from it.

**Nothing invents a coefficient.** Over thirty fits on genomes drawn with no runs of homozygosity at
all — five seeds at each of three budgets and two heterozygosities — every single one was refused.
That is the failure this measurement was most worried about, and it did not occur.

**Today the coefficient is not fitted at all**: `estimate-parameters --inbreeding` states one number
for the whole cohort and defaults it to zero, recorded as `supplied`. So the comparison that matters
is not 0.086 against 0.026; it is **0.086 against the 0.78 error a selfing panel gets from the
default**.

---

## 2. What was built, and why it is a bridge rather than a second estimator

There is one runs-of-homozygosity estimator in the tree, `generic::runs::fit_inbreeding`, and it
reads a table of *what each 100 kb window of the genome looked like*. The new module reshapes a
sample's census records into that table and changes nothing else, so a coefficient fitted from a
census and one fitted from the whole genome differ **only in which positions were looked at**. A
second implementation would have made them differ for two reasons at once and no comparison could
have said which it was seeing.

Its unit tests assert that over the same positions the two routes build the same cells — cell for
cell, over twelve positions landing in eight cells across two windows. **One real difference turned
up and it is not the bridge**: above the depth cap of 124 reads a position the whole-genome route
thins a position by a seeded hypergeometric draw while the census takes the proportional share, so
that a region-sharded walk writes byte-identical records. At 300 reads of which 150 disagree the
census records 62 and the draw returns 55. The census figure is the draw's expectation, so neither
is wrong — and the cap never fires on the tomato cohort at 2.5 to 28.6 reads a position.

---

## 3. What a budget buys: the accuracy sweep

**How to read it.** A census budget fixes how many positions a window holds, and the estimator
classifies a window on how many *heterozygous* positions are among them. On a tomato genome — 800 Mb
of analysable sequence, 8,004 windows of 100 kb — the two are one number:

| census positions | positions a window | heterozygotes a window at 0.91 per kb |
|---:|---:|---:|
| 2,000,000 (shipped, before and after) | 250 | 0.23 |
| **6,000,000** (measured here, not adopted) | 750 | 0.68 |
| 16,000,000 | 2,000 | 1.8 |
| 60,000,000 | 7,500 | 6.8 |
| every position | 100,000 | 91 |

Genomes are drawn with a known autozygous fraction at tomato's shape — 8,004 windows over twelve
contigs, three reads a position, a read misreading at 0.0033, homozygous-non-reference at 1.5 per
kilobase — and `F` is fitted with the shipped estimator. Five seeds a cell. `worst err` is the
largest `|fitted − realised|` over those five, which is the sound column: `realised F` varies by
seed, so a difference of means is not an error.

**A selfing line, 78% of the genome in runs of 3 Mb** — the tomato panel's own shape, whose cohort
fit put the median homozygote excess at 0.78:

| census budget | het a window | realised `F` | mean fitted | **worst error** | worst start spread | refused |
|---:|---:|---:|---:|---:|---:|---|
| 2M | 0.23 | 0.7797 | 0.7970 | **0.1054** | 0.302 | 1 of 5 |
| **6M** | **0.68** | **0.7756** | **0.7525** | **0.0856** | **0.0013** | **0 of 5** |
| 16M | 1.8 | 0.7832 | 0.7868 | **0.0256** | 0.0000 | 0 of 5 |
| 60M | 6.8 | 0.7729 | 0.7728 | **0.0039** | 0.0000 | 0 of 5 |
| every position | 91 | 0.7887 | 0.7887 | **0.0000** | 0.0000 | 0 of 5 |

**An outcrossing-to-moderate genome, 30% in runs of 3 Mb:**

| census budget | het a window | realised `F` | mean fitted | **worst error** | worst start spread | refused |
|---:|---:|---:|---:|---:|---:|---|
| 2M | 0.23 | 0.3295 | 0.3617 | **0.1108** | 0.0056 | 0 of 5 |
| **6M** | **0.68** | **0.3070** | **0.3261** | **0.0360** | **0.0001** | **0 of 5** |
| 16M | 1.8 | 0.2991 | 0.2986 | **0.0104** | 0.0000 | 0 of 5 |
| 60M | 6.8 | 0.3000 | 0.2989 | **0.0022** | 0.0000 | 0 of 5 |
| every position | 91 | 0.3241 | 0.3241 | **0.0001** | 0.0000 | 0 of 5 |

**The `worst start spread` column is how far apart the nine starting points the fit climbs from
landed**, largest coefficient minus smallest, over the worst of the five seeds. It is what the run
is accepted or refused on. Read it against §4's no-runs figures, which are two to three orders of
magnitude larger.

**Two things to take from these.** The error falls with the budget and does not plateau, so this is
a precision question rather than a broken-instrument one; and **six million is where the refusals
stop** at the cohort's median heterozygosity, which two million did not manage on the selfing arm.

---

## 4. Where it still refuses, and that is the right answer

Spec §8's fourth measurement asks for every figure split by the sample's own heterozygosity rather
than pooled, because at the cohort's floor a mis-fitted background shows up as a confident wrong
number rather than as scatter. The tomato panel's accessions run **0.32 to 4.83 heterozygotes per
kilobase, median 0.91**, and that spread moves the answer more than the budget does.

Selfing genomes at 78%, at the panel's two ends:

| census budget | het/kb | het a window | realised `F` | mean fitted | worst error | worst start spread | refused |
|---:|---:|---:|---:|---:|---:|---:|---|
| 2M | 0.32 | 0.08 | 0.7748 | — | — | 0.718 | **5 of 5** |
| 6M | 0.32 | 0.24 | 0.7898 | 0.8203 | 0.0707 | 0.194 | **3 of 5** |
| 16M | 0.32 | 0.64 | 0.7780 | 0.7503 | 0.0834 | 0.069 | **1 of 5** |
| 2M | 4.83 | 1.20 | 0.7732 | 0.7711 | 0.0246 | 0.0000 | 0 of 5 |
| 6M | 4.83 | 3.62 | 0.7774 | 0.7794 | 0.0053 | 0.0000 | 0 of 5 |
| 16M | 4.83 | 9.66 | 0.7653 | 0.7654 | 0.0009 | 0.0000 | 0 of 5 |

**At the panel's least heterozygous accessions six million refuses three samples in five, and
sixteen million still refuses one.** A refusal is not a gap: it is the estimator saying the starting
points it climbed from disagreed about the answer, which is what
`MAX_IDENTIFIED_START_SPREAD` exists to say. The design's own rule is that such a sample takes a
supplied coefficient rather than a fitted one, and the run names it.

**At the most heterozygous end, two million was already enough** — worst error 0.025. So the budget
is not a property of the cohort but of the sample, and what six million changes is how much of the
panel clears the bar.

### The failure that did not happen

`inbreeding_resolution_2026-08-09.md` §1 found the alarming case: on genomes with **no runs at all**
and about five heterozygotes a window, two fits in sixty returned `F ≈ 0.99` — a fabricated
coefficient, accepted, at twenty-four and a hundred and fifty times the reported noise floor. A thin
census sits squarely in that regime, so it was the thing to check.

**Thirty fits on genomes with no runs, across three budgets and two heterozygosities: thirty
refusals, and not one fabricated coefficient.** Read against that note, this is the
across-start-spread criterion — adopted after it, as spec §6.3 — doing the job it was adopted for. At
the whole-genome end the same genomes gave three refusals in five and a fitted 0.0017 to 0.0028 on
the two that were accepted, which is nothing and is correct.

**And the margin is not narrow.** The quantity the refusal is made on, the spread across the
starting points, separates the two populations by two to three orders of magnitude at every budget:

| worst start spread | genomes with runs | genomes with none |
|---|---:|---:|
| 2M | 0.0056 to 0.302 | **0.606 to 0.922** |
| 6M | 0.0000 to 0.194 | **0.492 to 0.989** |
| 16M | 0.0000 to 0.069 | **0.335 to 0.911** |
| every position | 0.0000 | 0.0061 |

The one overlap is the low-heterozygosity column of §4's table — 0.194 at 6M with real runs against
0.492 with none — and that is the same cell where three fits in five were refused. **So a thin census
does not produce a wrong coefficient; it produces no coefficient**, and the two failure modes are
not equally bad.

---

## 5. What it costs

Measured with `generate-psps` on four tomato accessions over the whole 8 Mb of
`benchmarks/tomato1/regions.bed`, at about three reads a position, before and after the budget
change. The command reports each psp's census as a share of the file, so this is read off the run
rather than computed.

| sample | census at 2M | census at 6M | increment |
|---|---:|---:|---:|
| SRS3394712 | 2,692,861 | 7,773,665 | +5,080,804 |
| SRS3394606 | 2,872,909 | 8,285,939 | +5,413,030 |
| SRS3394713 | 2,777,857 | 8,010,539 | +5,232,682 |
| SRS3394712_SRR7279484 | 2,477,162 | 7,127,706 | +4,650,544 |
| **mean** | **2,705,197** | **7,799,462** | **+5,094,265** |

**2.71 MB a sample to 7.80 MB — an increment of 4.86 MB, or 1.27 bytes a census position.** A byte
of that is the dense depth code every position carries; the remaining 0.27 is the sparse list of
what was not on the reference base.

**The whole psp grew by exactly its census and by nothing else** — 45.03 MB to 50.12 MB on the first
accession, and the same arithmetic on the other three — because the budget changes the trailer and
not the evidence.

**The size follows the budget and not the ground.** Six million positions cost the same 7.8 MB
whether the run covers 8 Mb or a whole genome, as long as there is more ground than budget; the
8 Mb benchmark is above it at both settings. So **on a 63-accession cohort the censuses go from about
171 MB to 491 MB in total**, and a whole-genome run costs the same.

### The fit is where the cost actually lands, and it is larger than the psp's

The psp figures above are the cheap half. `estimate-parameters` over the same four accessions,
peak resident memory from the kernel's own high-water mark (`scripts/peak_rss.sh`):

| | peak resident | wall clock |
|---|---:|---:|
| 2M | 185.8 MB | 944 s |
| **6M** | **516.9 MB** | **2,010 s** |
| | **2.78×** | **2.13×** |

**+331 MB and +18 minutes, on four samples.** That is about 83 bytes of peak memory a census
position — far more than the 1.27 bytes the position costs on disk — and **this report does not
attribute it**: the fit holds a per-position quantity for every kept position as well as the
samples' records, and which of those grew is not measured here.

**How it scales with the cohort is not measured either**, and the four-sample pair cannot say: the
fit's memory has a per-position part that does not grow with samples and a per-sample part that
does, and two budgets at one sample count cannot separate them. A 63-accession cohort should be
measured before this budget is relied on at scale.

### Walking the reads

Walking four accessions over the 8 Mb took **44.6 s at six million against 33.0 s at two million
for three of them**. The two runs are not a clean pair on sample count, so read this as *tens of
seconds either way on this ground* rather than as a ratio. The walk is not where the budget hurts.

### What tripling the census did to the parameters it was already fitting: almost nothing

The two runs' parameters files, on the same four accessions:

| fitted quantity | at 2M | at 6M | change |
|---|---:|---:|---:|
| genotype prior's reference concentration | 1.48750 | 1.48534 | **−0.15%** |
| genotype prior's alternative concentration total | 0.0041475 | 0.0041828 | **+0.85%** |
| the four read groups' error multipliers | 3.5246, 5.7364, 6.5206, 5.2912 | 3.5226, 5.7536, 6.5414, 5.2911 | **all within 0.3%** |
| fallback length-spectrum concentration | 66.379 | 64.595 | −2.7% |

**This is the strongest evidence in the report that two million was the right budget for everything
the census was previously asked for**, and it says the raise buys nothing for those parameters. What
it buys is the *possibility* of fitting the inbreeding coefficient, which §3 sizes — and both files
still record that coefficient as `supplied 0.0`, because the bridge of §2 is not wired into the
command yet.

---

## 6. What changing the budget does to psps already on disk

**Every census written under the old budget becomes stale**, because the position budget is one of
the settings a fit compares before it will pool two samples. That is not a silent hazard: running
the new binary against the old psps refuses, names all four samples, says which setting differs, and
prints the command that repairs them.

```
SRS3394606 (tmp/psp_2M/SRS3394606.psp) carries a census recorded under settings this run
does not use (the first that differs: generic target position count)
…
Rebuild the 4 samples whose census is named above, then run this fit again:
  regenerate-census --reference … --catalog … --psp tmp/psp_2M
```

`regenerate-census` rebuilds a trailer from the stored psp and **does not re-walk the reads**.

---

## 7. What this does not settle

- **The comparison is against a drawn genome, not against a real one.** There is no truth set for
  any real sample's autozygous fraction, which is why the whole-genome route's own accuracy was
  never established on real reads either. What a drawn genome can say is *does the estimator recover
  what it was given*, and what it cannot say is whether real runs look like drawn ones.
- **The sweep models the thinning as positions a window rather than by fitting real census
  records.** That rests on the bridge's unit test — same positions, same cells — which holds for
  every depth at or below the cap.
- **⚠ The budgets are not run on the *same* drawn genome, only on the same drawing process.**
  Drawing a window consumes a number of random values that depends on how many positions it holds,
  so the generator's state has diverged by the time the next window's run state is drawn, and each
  budget therefore gets its own run structure from a given seed. The `realised F` column shows it:
  0.3070 to 0.3295 across the budgets on the 30% arm, where a paired design would show one value.
  So each row is a fair estimate of *that budget's* error, and the rows are **not** a paired
  comparison — a difference of 0.02 between two adjacent rows is inside what the unpaired genomes
  contribute. The ordering across the whole range, 0.105 down to 0.0000, is far larger than that
  and is not at risk from it.
- **One heterozygosity per genome, and real accessions are not uniform.** A sample whose
  heterozygosity varies along the genome for reasons other than autozygosity is exactly what the
  two-state model is supposed to be fooled by, and nothing here draws one.
- **Nothing downstream was measured.** Whether a coefficient wrong by 0.086 changes a call is a
  separate question; the only sizing available is
  [`inbreeding_sensitivity_of_the_seed_2026-08-23.md`](../ng/reports/inbreeding_sensitivity_of_the_seed_2026-08-23.md),
  which measures errors of 0.05 to 0.10 moving the prior odds on a heterozygote by about 2% at 26
  to 63 individuals. An error of 0.086 is inside that range; the 0.78 the default carries is not,
  and that report says nothing about an error that size.
- **The bridge is not wired into `estimate-parameters` yet.** The command still declares the
  coefficient. Making it fit one is the next step and is where the per-sample refusals above have to
  be surfaced.
- **Sixteen million was not costed.** At the measured 1.27 bytes a position it would be about
  20.3 MB a sample, but that is arithmetic rather than a measurement.
