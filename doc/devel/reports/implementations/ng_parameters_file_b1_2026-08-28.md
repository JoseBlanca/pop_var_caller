# ng parameters file — B1: `RunParameters` → the file shape

**Date:** 2026-08-28
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone B, step B1
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §3
**Code:** [src/ng/calling/parameters_file/from_run_parameters.rs](../../../../src/ng/calling/parameters_file/from_run_parameters.rs), with read-only accessors added to [src/ng/calling/run_parameters.rs](../../../../src/ng/calling/run_parameters.rs) and [src/ng/parameter_estimation/joint/stratum_fits.rs](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs)

---

## 1. Plan

Map every number a run scored its reads under onto the row of the file that carries it, with the
read-group and sample axes carrying their names as well as their indices. No TOML: this step
produces a `ParametersFile` value and stops, which is what lets the projection be tested on the
shape rather than on a rendering of it.

The entry point is `ParametersFile::of_run`, in a new submodule beside the shape.

## 2. The projection reads more than `RunParameters`, and that was A2's finding

Assembly drops two things spec §3 asks for by name, so a projection written from the assembled
parameters alone would produce a plausible file that is wrong about both:

- **the base-quality calibration's evidence count** (§3.3). `ReadGroupCalibration` is a multiplier
  and a warrant; the count of reads is on the `Estimate<ErrorRate>` that `RunParameters::assemble`
  reads and does not store.
- **every sample's inbreeding warrant and count** (§3.5). The pre-pass fits an
  `Estimate<InbreedingF>` per sample; the seam takes a bare `Vec<InbreedingF>`, so a file written
  from it would mark every fitted coefficient as handed over.

Both estimate sets are therefore arguments, beside the run's read-group table — which is where the
sample names, the `@RG ID`s and the library names live, none of which the parameters carry — and
beside the reference digest and the census identity, which are what the file is bound to (§3.1)
and are not derivable from any number.

**Six arguments, no bundle type.** `RunParameters::assemble` takes nine and its own documentation
says why: the list is the point, and a struct naming the same six things is a second place for
them to go out of step with the run.

## 3. Assumptions and choices, none of which the spec settles

- **A `defaulted` value carries no evidence count.** Read group 1 of the fixture fitted a rate of
  zero over 4,242 reads; `ReadGroupCalibration::from_fitted_rate` refuses it and the calibration
  becomes the honest defaulted one, multiplier 1.0. Writing 4,242 beside that 1.0 would say the
  multiplier rests on those reads, and it rests on nothing. The two other defaulted numbers in the
  file — the tract ladder's stated concentration and the outlier weight — already carry no count,
  so this is one rule everywhere rather than a rule about calibration.
- **A slippage number that came off a curve must carry the curve, or the projection panics.** The
  file's `LevelSmoothing` and `ShareSmoothing` hold the curve and the reach *on the variant*, by
  A1's design, so there is no legal shape for "the curve supplied this and I cannot say which".
  Falling back to *the stratum's own fit* would write a false warrant; a panic naming the stratum
  and the slippage group is the alternative. **Unreachable from the fit** as it stands — the
  smoothing pass sets source and curve together — and every field of the upstream provenance is
  public, which is why it is checked rather than assumed.
- **`LevelSource::Cell` ignores any curve beside it.** That source means the stratum's period had
  no curve at all (`blend_level` returns it only for a fitted cell with no curve), so there is
  nothing to lose by not reading one.
- **A read group the slippage declaration does not name gets no row** in
  `slippage_group_by_read_group`, which is what the fit itself says about it: `StratumFits::at`
  answers `UnknownReadGroup` and no slippage number is ever looked up under it.
- **The outlier weight is written `defaulted` from the compiled-in constant**, because
  `RunParameters` holds no other value — see §8 for the gap that leaves.
- **Every mismatch check is a release panic**, mirroring `RunParameters::report`'s pair: two
  tables minted from different inputs and joined positionally write one library's numbers under
  another's name, which looks like an answer. **The review disagreed** — this runs after the last
  locus, so a panic discards a cohort's calling work where `assemble`'s equivalent checks discard
  a startup. The panics stay, because each condition is a wiring bug in whoever assembles the six
  arguments rather than a state any input data reaches, and the doc now argues that case rather
  than pointing at `assemble`'s. **The decision needs a call site and there is none yet**; step F1
  decides the order of the two writes and therefore whether a failed projection can cost a VCF.

## 4. A gap in the shape that this step found: the seed had three rungs and needs four

`SeedRegime` has **four** states and the file's `SeedRung` was built with three. The missing one
is `ZeroDiversity` — the run's measured heterozygosity was exactly zero, a cohort with no
variation at all. Its pair of concentrations is a legal pair that says nothing about how it was
arrived at, so the rung is the only place that state is recorded; a run there would have had no
rung to be written under.

It was added, with the fourth spelling pinned in
`every_enum_variant_spells_as_the_file_says`. **This is the drift `PROJECT_STATUS.md` flagged as
unguarded after A1**, and the guard is now structural rather than a checklist: every one of the
eight conversions from a pre-pass enum to its file spelling is an exhaustive `match`, so a variant
added upstream fails to compile here.

**One number in A1's prose was wrong and is corrected in the same commit.** The doc on
`every_enum_variant_spells_as_the_file_says` said the golden file exercises "eleven of these
twenty-one" unit variants. Counted by grepping the golden file for each asserted spelling: it is
**fourteen of twenty-two** (twenty-one before this step's fourth rung). The eight it never writes
are the seed's other three rungs, three of the share curve's four, and two of the share shape's
three.

## 5. Changes made

| file | change |
|---|---|
| `parameters_file/from_run_parameters.rs` | new — `ParametersFile::of_run`, its five section builders, and the eight enum conversions |
| `parameters_file/mod.rs` | `SeedRung::ZeroDiversity`; `mod from_run_parameters;`; the corrected count in one test's doc |
| `calling/run_parameters.rs` | seven read-only accessors: `ploidy`, `calibration_by_read_group`, `contamination_by_read_group`, `inbreeding_coefficient_by_sample`, `prior_seed`, `ssr_slippage_fits`, `ssr_substitution_rate`. They give out what is stored and compute nothing |
| `parameter_estimation/joint/stratum_fits.rs` | three iterators — `each_stratum_and_group_with_numbers`, `fitted_length_spectrum_of_each_stratum`, `pooled_length_spectrum_of_each_period` — and the shared `fitted_at` helper they use with `at`. `StratumFits::at` answers one lookup and cannot enumerate, and writing the run's parameters down means walking every cell. **The names say the walk is not dense**: two of the three skip the cells with no numbers, and a caller that zipped them against a per-stratum table would misalign it |

## 6. Tests

**After the review's fixes: 31 tests in the projection module and 3 more in `stratum_fits.rs`.**
The first draft had 16, of which the review proved one could not fail and six could not see a
mutation; the counts here are the fixed tree's, and the fix report lists the six mutations re-run
against it.

The fixture is a run assembled through `RunParameters::assemble` rather than built by hand: three
read groups over two samples — two lanes of one plant and one of another, so the two per-axis
lengths differ — one contaminated read group, one that identified nothing and one measured and
found clean, one fitted stratum, two derived from their period's curves and one refused, and two
substitution rates at different read groups, periods and **ploidies**.

The tests that carry the step:

- `the_calibration_carries_the_count_the_assembled_parameters_dropped` and
  `each_samples_inbreeding_carries_the_warrant_the_seam_dropped` — the two things a projection
  from `RunParameters` alone cannot supply.
- `a_defaulted_number_carries_no_evidence_count` — discriminating because the refused rate behind
  that multiplier of one has 4,242 reads on it, so the absence cannot come from a zero count. It
  checks the same rule on a defaulted substitution rate, so the rule is the file's and not the
  calibration's. `an_evidence_count_of_zero_is_written_and_is_not_absence` holds the other
  direction: a count of zero is a count.
- `an_unmeasured_read_group_writes_no_measurement_and_a_clean_one_writes_a_zero`, plus
  `an_uncontaminated_run_writes_no_contamination_section` and
  `a_pair_with_no_reads_and_a_stratum_with_no_fit_get_no_row` — four of spec §5's five states, on
  the way out. (The fifth, a defaulted multiplier of 1.0, is the defaulted-count test's fixture.)
- `every_pre_pass_word_maps_to_its_own_word_in_the_file` — the other half of the drift guard.
  Exhaustiveness catches a variant *added* upstream; nothing but this table catches one *crossed*
  with another, and `Flat => Sloping` compiles and spells correctly.
- `of_run_writes_a_single_sample_single_read_group_run` — the bottom of the committed input range,
  and the only shape in which the read-group and sample axes are the same length.

**One assertion is against a bound rather than an equality, and the bound is the accumulator's.**
The fitted multiplier is a fitted rate of 0.004 over reads averaging 0.008, which is a half — and
comes out as **0.4999999998066437**, because the denominator is summed in fixed point so that
shards merged in different orders give the same number. That is **3.9 in 10¹⁰ below a half**,
against the 4.8 in 10⁷ `ReadGroupCalibration`'s own documentation states for the mean it is built
from. The test asserts against that stated bound and says where it comes from.

## 7. Validation

Run in the dev container, `./scripts/dev.sh`:

| command | result |
|---|---|
| `cargo fmt --check` | clean, exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean, exit 0 |
| `cargo test --lib ng::calling::parameters_file` | **45 passed, 0 failed, 1 ignored** |
| `cargo test --lib ng::parameter_estimation::joint::stratum_fits` | **30 passed, 0 failed** |
| `cargo test --lib` | **4,968 passed, 0 failed, 12 ignored** |

`cargo test --all-targets --all-features` is **not** the gate here: it exits 101 on a pre-existing
panic in `benches/psp_writer_perf.rs:386`, verified on clean `main` and recorded in
`PROJECT_STATUS.md` as a standing project-wide item.

## 8. Follow-ups this step surfaced

- **An edited outlier weight has nowhere in memory to live.** Spec §3.8 says a person editing that
  number is the point of writing it down, and its `supplied` state exists for exactly that — but
  the weight reaches calling from `likelihood/ssr.rs` rather than from `RunParameters`, so a run
  assembled from a supplied file would ignore the edit and this projection would write the
  compiled-in constant back. **Step C2's**, where the reader's projection runs.
- **The projection can still produce a contamination table in which no row has a measurement**, if
  every read group's estimate came back `Estimated` with zero counts. That is one of the two
  shapes `the_shape_accepts_two_things_step_c2_must_refuse` pins, and C2 refuses it. The review
  probed for it and **could not reach it from the estimator** — `fit_contamination_over` refuses
  below 100 markers — only by hand. It is left faithful rather than collapsed at the writer: in
  that state the run took the *mixture* path with every fraction zero, and absence asks a reader
  for the plain formula instead, so collapsing would make the file say something the run did not
  do.
- **The design-fidelity pass found one spec sentence the code contradicts**, offered for the spec
  rather than the code: §3.6 says an alternative concentration total of exactly zero "is not
  floored on the way in or out", and the seed builder floors it before `RunParameters` ever holds
  it. `SeedRung::ZeroDiversity` — added by this step — is what carries that state instead.
