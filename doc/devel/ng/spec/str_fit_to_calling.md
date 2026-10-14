# What the repeat-tract fit hands to calling

**Status:** draft for the owner's approval, 2026-10-14. Companion of
[`str_penalised_fit.md`](str_penalised_fit.md), which settles the fit; this settles what reaches
calling from it: the parameters file's repeat-tract section, and the fallback calling uses when a
class has no between-locus share. Build the two together.

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
SNP/indel sections of the file. It does not compute the period-pooled rung (§3.3).

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
| `fallback_length_spectrum_concentration` | `fallback_between_locus_share` | F of the fallback; warrant kept |

`FORMAT_VERSION` becomes 3 (`parameters_file/mod.rs:192`). Validation requires F strictly between 0
and 1 (an F of 0 is an infinite concentration, 1 a zero one), in place of a positive finite
concentration.

**Convert at the file's edge only.** The fit and calling keep working in the concentration; one pair
of functions converts, and nothing else touches the formula.

### 2.2 Files written before version 3 still read

The project reads older files (`mod.rs:187-191`: a version-1 file reads as a version-2 one with the
new keys absent). **Version 1 and 2 files are migrated on the parsed TOML before they are
deserialised**, in one function in `from_toml.rs`: the keys of §2.1 renamed, each `concentration = c`
replaced by `between_locus_share = 1/(1+c)`, the fallback's value converted and its warrant kept. A
non-positive or non-finite old concentration is not converted into a valid share; it fails validation
as before. **This is the half of the renaming the branch lacks**: today it refuses every older file,
which the kimura session had to work around by hand (`research/str_fit/kimura_shapes_span9_2026-10-14.md`).

### 2.3 Where a stutter number came from

Each stratum's slip rate, shorter share and fall-off carry their origin (`LevelOrigin`,
`parameters_file/mod.rs:1125`): main's fit writes the stratum's own fit, its period's parametric curve
(`SlippageCurve`), or a blend of the two (`LevelSmoothing`, `:1190`; `LevelSource`,
`slippage_curve.rs:519-531`). The new fit produces none of these: every number is a point on a
penalised curve, fitted with the class's own reads.

**Version 3 adds one variant, `PenalisedCurve`**, carrying: the curve's stiffness strength; whether
this repeat count lies inside the range of classes with reads (outside, the value is held at the
nearest end); whether a bend was eased at this repeat count (a kink the data insisted on); and the
expected slipped reads behind the class. `own_fit_standard_error` holds the fit's error
(`str_penalised_fit.md` §2.6). The old variants stay readable for older files. *A variant rather than
reusing `ThisStratum`: a value pulled toward its neighbours is not the stratum's own fit, and the file
says which.*

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
measured:* a small cohort where F is absent because its estimate failed, which the new fit does not
report as absent; and an outbred single sample.

**Placed above the pooled step** because a class's own shares were measured better than flat, and the
pooled step is not produced in production (§3.3), so its order against the new step is untested.

**What changes in code:** the class's concentration in calling's tables becomes optional; the step is
a new `AlleleLengthVariationRung` value (`stratum_fits.rs:359-366`), recorded per locus as the rung
is today (`summarise_condition.rs:2396-2447`). The rung is not written to the VCF today
(`vcf/assemble.rs:103`), so no output format changes.

### 3.3 The pooled step stays as it is

It is filled only from a parameters file that carries `allele_length_variation_by_period`; no
production run computes it (`fit_period_allele_length_variations`, `ssr_fit.rs:3390`, has only test
callers). Unchanged here.

### 3.4 The stated F

The median of the fitted classes' F, converted from the median concentration (the conversion is
monotone, so the median carries over), or F = 0.5 (concentration 1) where no class has one. At one
sample no class has one, so every class uses 0.5.

---

## 4. What will bite

- **The two specs land together.** The new fit returns F absent at one sample; calling panics on it
  until §3.2 lands.
- **The golden files** (`parameters_file/testdata/every_shape_as_written.toml`, `every_shape.toml`)
  are regenerated by an ignored test (`to_toml.rs`, `regenerate_the_written_golden_file`); check the
  diff shows only the intended keys. **`version_1_as_written.toml` must keep the old keys**: the
  branch's mechanical rename wrongly renamed them (`aa42d629`); restore it from the commit before.
- **The oracle** re-records the fitted file and the calls made with it; the default-parameter outputs
  should not move.

## 5. How we know it works

1. A version-1 and a version-2 file read back with every share's F equal to 1/(1+c) of the old
   concentration, and every other number unchanged.
2. A version-3 file round-trips through write and read unchanged, F absent included.
3. Validation refuses F of 0 or 1.
4. A class with shares and no F is seeded from its own shares at the stated concentration; one with
   neither falls to flat (tests in `seed_ssr.rs` and `stratum_fits.rs`).
5. On the bench at one sample, calling's rung for every class with reads is the new step, and the
   held-out score matches the measurement of §3.2.

## 6. Open

1. **Should the file carry errors for the allele length variation and F**, as it does for the stutter
   numbers? The fit computes them (`str_penalised_fit.md` §2.6). *Leaning:* F's error yes, the shares'
   not yet, because those for rare lengths are not honest (`str_penalised_fit.md` §5.2). Confirm
   before code.
