# What the repeat-tract fit hands to calling

**Status:** approved by the owner, 2026-10-14. Companion of
[`str_penalised_fit.md`](str_penalised_fit.md), which settles the fit; this settles what reaches
calling from it: the parameters file's repeat-tract section, and the fallback calling uses when a
class has no between-locus share. Build the two together.

**Amended 2026-10-14 after an adversarial review.** Corrections of fact are made in place; **the
lead's decisions are marked *(lead, 2026-10-14)* and listed in §6 for the owner to confirm.**

**Amended 2026-10-15 by the owner:** the stated F's median counts only classes with at least 8 loci
with reads (§3.5).

---

## 1. What this is

The fit writes its numbers into a parameters file (TOML, `calling/parameters_file/`), and calling
reads them back into its tables (`StratumFits`, `parameter_estimation/joint/stratum_fits.rs`). Three
things in that hand-over change:

1. **Names and the number reported for spread.** The owner's renaming of 2026-10-12: "length
   spectrum" becomes **allele length variation**, and the Dirichlet concentration the maths uses is
   reported as the **between-locus share** F = 1 / (1 + concentration), the Balding–Nichols F_ST with
   a class's loci in place of populations (0: every locus carries the class's whole variation; 1:
   each locus is fixed at one length).
2. **F may be absent.** The new fit reports no F where nothing in the data can place it, at one
   sample (`str_penalised_fit.md` §2.5). Today calling cannot hold that.
3. **Where a stutter number came from.** The new fit's curves are not main's curves, and the file's
   provenance fields assume main's.

**Goals.** Calling gets every number the new fit produces, in the owner's names; a class with its own
allele length variation but no F keeps its own variation; a file written by an earlier build still
reads.

**Non-goals.** It does not change how calling scores reads, which candidates it considers, or the
SNP/indel sections of the file. It does not compute the period-pooled rung (§3.4).

---

## 2. The parameters file, format version 3

### 2.1 The keys

The renaming of the two tables is already on the branch (`aa42d629`); the rest is new.

| version 2 | version 3 | value |
|---|---|---|
| `length_spectrum_by_stratum` | `allele_length_variation_by_stratum` | same rows |
| `length_spectrum_by_period` | `allele_length_variation_by_period` | same rows |
| row `shares_by_repeat_offset` | row `share_of_chromosomes_by_repeat_offset` | same numbers |
| row `concentration` | row `between_locus_share` | F = 1/(1+c); **absent where F is** |
| — | row `between_locus_share_standard_error` | F's error; absent where F is |
| `fallback_length_spectrum_concentration` | `fallback_between_locus_share` | F of the fallback; warrant kept |

`FORMAT_VERSION` becomes 3 (`parameters_file/mod.rs:192`). Validation requires F strictly between 0
and 1 (an F of 0 is an infinite concentration, 1 a zero one), in place of a positive finite
concentration, **and requires F's concentration to be finite**: an F below about 5.6×10⁻³⁰⁹ passes
"strictly between 0 and 1" but converts to an infinite concentration, which calling asserts against
(`AlleleLengthVariation::fitted`, `stratum_fits.rs:257-262`).

**Convert at the file's edge only.** The fit and calling keep working in the concentration; one set of
functions converts, and nothing else touches the formula.

**F's error crosses on the concentration's scale** *(lead, 2026-10-14)*. The fit and calling carry the
concentration's error (`StratumErrors::concentration`, `ssr_fit.rs:611-621`); the file carries F's,
converted at the edge, beside the two conversions, by the delta method: se(F) = se(c) / (1 + c)².

### 2.2 Files written before version 3 still read

The project reads older files (`mod.rs:185-191`: a version-1 file reads as a version-2 one with the
new keys absent). **Version 1 and 2 files are migrated on the parsed TOML before they are
deserialised**, in one function in `from_toml.rs`: the keys of §2.1 renamed, each `concentration = c`
replaced by `between_locus_share = 1/(1+c)`, the fallback's value converted and its warrant kept.
**This is the half of the renaming the branch lacks**: today it refuses every older file, which the
kimura session had to work around by hand (`research/str_fit/kimura_shapes_span9_2026-10-14.md`).

What the migration does, exactly *(lead, 2026-10-14)*:

- **Both spellings of version 2 are accepted**: main's tables (`length_spectrum_by_stratum`,
  `length_spectrum_by_period`) and this branch's since `aa42d629` (`allele_length_variation_by_*`,
  still with `concentration` in their rows, as in the oracle's `fitted.parameters.toml`). A row's
  `concentration` is converted whichever table name holds it.
- **The migrated file says version 3.** `format_version` is set to 3 by the migration, so a file read
  and written back carries version-3 keys under a version-3 stamp.
- **A version-3 file is deserialised from its own text**, so its errors name a line as today
  (`from_toml.rs:90-101`). **A version-1 or version-2 file is deserialised from the migrated TOML,
  which has no positions: its errors name the key path, not a line.**
- **An old concentration that cannot become a valid F is refused**, not converted into one: a
  non-positive or non-finite one, and one so small (below about 1.1×10⁻¹⁶) that F rounds to 1.

### 2.3 Where a stutter number came from

Each stratum's slip rate, shorter share and fall-off carry their origin (`LevelOrigin`,
`parameters_file/mod.rs:1134`): main's fit writes the stratum's own fit, its period's parametric curve
(`SlippageCurve`), or a blend of the two (`LevelSmoothing`, `:1199`; `LevelSource`,
`slippage_curve.rs:519-531`). The new fit produces none of these: every number is a point on a
penalised curve, fitted with the class's own reads.

**Version 3 adds one variant, `PenalisedCurve`**, carrying: the curve's stiffness strength; whether
this repeat count lies inside the range of classes with reads (outside, the value is held at the
nearest end); whether a bend was eased at this repeat count (a kink the data insisted on); and the
expected slipped reads behind the class. The old variants stay readable for older files. *A variant
rather than reusing `ThisStratum`: a value pulled toward its neighbours is not the stratum's own fit,
and the file says which.*

**`own_fit_standard_error` holds the fit's error beside a `PenalisedCurve` number**
(`str_penalised_fit.md` §2.6). Two places decide today whether such an error is kept, and both
accept only the stratum's own fit; both must accept `PenalisedCurve` as well:

- `OwnFitStandardErrors::of_the_numbers_emitted` (`stratum_fits.rs:107-126`) keeps an error only beside
  `LevelSource::Cell` / `ShareSource::Stratum`, so every error of the new fit would be dropped before
  it reaches the file, with nothing failing;
- the validator's `an_own_fit_error` (`parameters_file/validate.rs:1184-1200`, called at `:810-817`
  and `:851-865`) refuses an error beside any smoothing but `ThisStratum`, so a file carrying them
  would be written and then refused by every calling run.

**Calling's warrant for a `PenalisedCurve` number is `fitted_here`** *(lead, 2026-10-14)*: every class
the new fit reports has reads of its own (`str_penalised_fit.md` §2.2). Today the warrant is read off
the curve's weight (`warrant_of`, `calling/inference/repeat_tract_parameters.rs:1154-1165`); a
penalised point has no such weight.

---

## 3. The fallback ladder in calling

### 3.1 Today

For each repeat tract, calling looks up its class's allele length variation and the concentration it
is held with (`allele_length_variation_at`, `stratum_fits.rs:781-805`):

1. the class's own fit, shares and concentration together (`FittedAlleleLengthVariation`,
   `:379-385`);
2. else its motif length's pooled fit;
3. else a flat distribution at the stated concentration: the median of the fitted classes'
   concentrations, or 1 (`STATED_FLAT_CONCENTRATION`, `:409`) where none was fitted (`:1100-1112`).

**A class has shares only together with a concentration.** A class whose F is absent falls to the
flat bottom step, and a non-positive concentration panics (`:1051-1055`).

### 3.2 The new step — decided by the owner, 2026-10-14

**A class with its own allele length variation but no F keeps its own variation, held at the stated
concentration**, in place of falling to the flat step. The ladder becomes:

1. own variation, own F;
2. **own variation, stated F** (new);
3. the motif length's pooled fit;
4. flat, stated F.

**Measured** on the bench at one sample (`research/str_fit/track1_counting_fit_errors_2026-10-12.md`
§13): held-out prediction 0.11–0.16 log units a locus below the truth with own variation kept, against
0.84–1.68 with the flat step. **The stated value barely matters there**: 0.5, 1, 2 and 5 moved the
result by at most 0.02 a locus, because with one sample the concentration only changes how often a
sample's two copies match, and the bench's samples were highly inbred (excess 0.85–0.99). *Not
measured:* a small cohort where F is absent because its estimate failed; and an outbred single
sample.

**Placed above the pooled step** because a class's own shares were measured better than flat, and the
pooled step is not produced in production (§3.4), so its order against the new step is untested.

**What changes in code:** the class's concentration in calling's tables becomes optional; the step is
a new `AlleleLengthVariationRung` value (`stratum_fits.rs:359-366`), recorded per locus as the rung
is today (`summarise_condition.rs:2453-2526`), reached through a new `FittedFrom` value
(`stratum_fits.rs:190-195`). **The rung's position in its enum is its rank**: the enum derives `Ord`
(`:358`) and a bundle reports the largest rung of its runs as its weakest
(`BundlePrior::rung`, `calling/genotype_prior/seed_ssr.rs:374-380`), so the new rung sits between the
own fit and the pooled step, in ladder order. The rung is not written to the VCF today
(`vcf/assemble.rs:107`), so no output format changes.

### 3.3 Repeat-tract bundles follow the same rule

A bundle (several runs called as one locus, `calling/genotype_prior/seed_ssr.rs`, merged from main
on 2026-10-14) is seeded run by run from each run's class. A run whose length falls past the fit's
window reads its edge class, shared among the distinct lengths that run takes there, exactly as a
single segment does (`run_lengths_at_edge`). Main gave such a run a weight of 0; the merge changed it
so that a one-run bundle is still seeded exactly as its segment, which main's own test requires. The
new step of §3.2 applies to each run's class the same way.

### 3.4 The pooled step: readable, no longer produced

It is filled only from a parameters file that carries `allele_length_variation_by_period`
(`StratumFits::of_gathered_rows`), and **it stays readable**. Its only producer,
`fit_period_allele_length_variations`, had only test callers and ran the climb; **it is deleted with
the climb, and the tests that used it build their period tables by hand** *(lead, 2026-10-14)*.
*Done 2026-10-16:* deleted with the climb, together with the builder that handed its output to
calling (`with_period_allele_length_variations`) and its type (`PeriodAlleleLengthVariation`); the
tests now pass their hand-built tables through the file's own door, `of_gathered_rows`.

### 3.5 The stated F

**The median of the fitted classes' F, taken in F** *(lead, 2026-10-14)*: at an even count, the mean of
the middle two F values, not F of the mean of the middle two concentrations, which differs. F = 1/11
(concentration 10, `STATED_FLAT_CONCENTRATION`) where no class has an F. Calling holds it as its
concentration.

**Which classes the median counts** *(owner, 2026-10-15)*: **only classes with an F and at least
`LOCI_FOR_AN_ERROR` (8) loci with reads.** If no class has both, the median is taken over every class
with an F, as before 2026-10-15. Where no class has an F, F = 1/11 as above.

Why: on the tomato2 oracle (4 accessions, 20 regions, about 5 reads a sample), the milestone C review
(`doc/devel/reports/reviews/str_penalised_fit_milestone_C_2026-10-15.md`) found that the median
moved from 0.526 under the climb to 0.723 under the penalised fit, and that about half of that move
came from which classes were counted, not from the F values themselves. The climb, which refused
classes below 8 loci, wrote F for 15 classes. The penalised fit, which fits every class with reads,
writes it for 46, and 31 of those have fewer than 8 loci. A thin class's F is mostly its motif
length's curve, and for some (three bases at 30 repeats: one locus, one sample) it is the curve
extended past the last class with a pair of samples; yet each such class got one vote, the same as
the one-base class at 8 repeats with 1,214 loci. Over the 15 classes with 8 or more loci the median
is 0.615.

The threshold is the same constant below which a class's F, slip rate and allele-length shares carry
no error (`str_penalised_fit.md` §2.6), so a change to `LOCI_FOR_AN_ERROR` moves it too. *2026-10-16,
owner:* after plan step D1 the shorter share and the fall-off got their own, lower line of 3 loci
(`LOCI_FOR_A_DIRECTION_ERROR`); `LOCI_FOR_AN_ERROR` stays 8, and so does this median's threshold.

**Its warrant says whether it was measured** *(lead, 2026-10-14; restated 2026-10-15 for the owner's
rule above)*: `fitted_here` whenever the median is taken from any class's F, whether over the classes
with 8 or more loci or, failing those, over every class with an F; `defaulted` only when no class has
an F and the stated 1/11 is used. So the test is still "does any class have an F", and the 2026-10-15
rule changes the value, not the warrant. A median taken only from thin classes is written as
`fitted_here` like any other; the file does not say how many loci stood behind it.

**The stated value was F = 0.5 (concentration 1) until 2026-10-17, and the owner moved it to 1/11.** At
one chromosome every homozygote is more probable than the heterozygote — 1.5 times at two candidate
lengths, `(K + 1)/2` at K — which at one outbred sample tipped near ties to the wrong homozygote
(a colleague's GIAB chr22 benchmark: 167 records at 5x, 56 at 10x). Measured on the GIAB per-sample
set, one sample at a time with `--defaults`, the move makes 17 repeat-tract genotypes right and 1 wrong
over 5x–30x and moves nothing at 300x; a fitted cohort never reads the value and is unchanged
(`../../reports/ng_str_single_sample_prior_2026-10-17.md`). *Not settled:* a one-sample run with its
own fit, where this value holds the class's reference-heavy shares, loses 4 genotypes net at 5x (same
report, §3.4).

Before step A5 of the plan the warrant was set by whether any class had an allele length variation,
so at one sample, where every class with reads has shares and none has F, the compiled-in 0.5 would
have been written as `fitted_here`, and the validator, which checks a `defaulted` value only against
the constant (`validate.rs:750-768`), could not catch it. It now asks whether any class has an F
(`StratumFits::over`, `stratum_fits.rs:657-674`).

**One value for every motif length.** This is today's rule (one median over every stratum,
`median_concentration`, `stratum_fits.rs:1250`), inherited and not measured against a median per motif
length; soft.

---

## 4. What will bite

- **The two specs land together.** The new fit returns F absent at one sample; calling panics on it
  until §3.2 lands.
- **The golden files** are regenerated by two ignored tests: `regenerate_the_golden_file`
  (`parameters_file/mod.rs:2037`, for `every_shape.toml`) and `regenerate_the_written_golden_file`
  (`to_toml.rs:1559`, for `every_shape_as_written.toml`); check each diff shows only the intended keys.
  **`version_1_as_written.toml` must keep the old keys**: the branch's mechanical rename wrongly
  renamed them (`aa42d629`); restore it from the commit before.
- **Tests that pin today's versions flip.** `a_version_this_build_does_not_read_is_refused`
  (`validate.rs:1733-1739`) asserts version 3 is refused; `from_toml.rs:367` asserts a version-1 file
  reads back as version 1, where the migration now stamps 3 (§2.2).
- **A benchmark script writes the old key**: `benchmarks/ssr_hg002/src/sweep_tract_parameters.sh:136`
  rewrites `fallback_length_spectrum_concentration` with a concentration; on a version-3 file its own
  check fails. It is changed to write `fallback_between_locus_share` with an F.
- **The oracle**: the VCFs called with default parameters should not move. The parameters files'
  checksums do: `from_alignments.parameters.toml` (format version, the renamed keys, the fallback as
  F) and `fitted.parameters.toml` (the same, plus every row's F). **The VCF called with the fitted
  file can move by rounding**, because the concentration read back is 1/F − 1 of a rounded F, not the
  fitted value bit for bit; such a move is traced to that round trip before the checksum is
  re-recorded.

## 5. How we know it works

1. A version-1 file and both spellings of a version-2 file read back as version 3, with every share's
   F equal to 1/(1+c) of the old concentration, and every other number unchanged.
2. A version-3 file round-trips through write and read unchanged, F absent included.
3. Validation refuses F of 0 or 1, and an F whose concentration is not finite; the migration refuses an
   old concentration that would give F = 1.
4. A class with shares and no F is seeded from its own shares at the stated concentration; one with
   neither falls to flat (tests in `seed_ssr.rs` and `stratum_fits.rs`). The new rung ranks between
   the own fit and the pooled step.
5. A `StratumFits` built from the new fit's outcomes at one sample puts every class with reads on the
   new step, and its stated F is `defaulted` (§3.5).
6. A fitted file round-trips with its stutter numbers' errors and F's error present (§2.3).
7. The stated F's median skips classes with fewer than 8 loci with reads when any class has 8 or
   more: a cohort whose thin classes would move the median gives the median of the thick classes
   alone, `fitted_here`; a cohort with no class of 8 loci gives the median of every class with an F,
   `fitted_here` (§3.5).

## 6. Decided

- **2026-10-14, owner: the file carries F's standard error** (`between_locus_share_standard_error`,
  absent where F is), **and not yet the allele lengths' errors**, because those for rare lengths are
  not honest (`str_penalised_fit.md` §7).

### Decided by the lead in response to the review, 2026-10-14 — for the owner to confirm

- F's error is carried on the concentration's scale and converted at the file's edge (§2.1).
- The migration accepts both version-2 spellings, stamps version 3, and reports an old file's errors
  by key path (§2.2).
- A `PenalisedCurve` number keeps its error and is `fitted_here` in calling (§2.3).
- The pooled step stays readable and loses its producer with the climb (§3.4).
- The stated F is the median taken in F, and its warrant is `defaulted` where no class has an F (§3.5).

### Decided by the owner, 2026-10-15

- **The stated F's median counts only classes with at least `LOCI_FOR_AN_ERROR` (8) loci with
  reads**; if none has 8, every class with an F, as before (§3.5). Its warrant stays `fitted_here`
  whenever any class has an F *(lead's reading of the rule)*. Built in the plan's step C9.
