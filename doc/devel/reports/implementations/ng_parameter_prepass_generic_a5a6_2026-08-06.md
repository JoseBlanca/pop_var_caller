# ng step 4, generic path — A5+A6: the output types, the error, and the fit floors

**Date:** 2026-08-06
**Plan:** [parameter_prepass_generic.md](../../ng/impl_plan/parameter_prepass_generic.md), Milestone A steps A5 and A6
**Design authority:** [arch](../../ng/arch/parameter_prepass_generic.md) §2.4, §5.2, §5.3, §5.4 · [spec](../../ng/spec/parameter_prepass_generic.md) §5, §6, §7

---

## 1. Plan

The last two steps of Milestone A, run as one loop and bundled deliberately: A5 declares
what the path emits and A6 declares what it emits *instead* when there is not enough
data. They are two halves of one surface, and neither is more than declarations.

- **A5** — `Provenance`, `Estimate<T>`, `SampleRates` with its two diploid accessors,
  `GenericSampleParameters`, `FitTermination`, `ScanResult<P>`, `CoupledFit`,
  `RunsModelStarts`, `RunsModelFit`, `StartOutcome`.
- **A6** — `ParameterEstimationError`'s four variants, with `MIN_SITES_TO_FIT`,
  `MIN_WINDOWS_TO_FIT_INBREEDING`, `DEFAULT_ERROR_RATE` and
  `MAX_COUPLED_FIT_ITERATIONS`.

## 2. Assumptions

Four choices the design left open. Each is a placement or a signature, none is a
decision about what the estimator does.

1. **Where each type lives**, which the architecture gives as one flat list. Split by
   who reads it, following the A1–A3 review's ruling that `parameter_estimation/mod.rs`
   is the level the STR sub-unit will share:
   - `mod.rs` — `Provenance`, `Estimate<T>`, `ParameterEstimationError`. Every parameter
     step 4 emits carries a provenance and an observation count, on either path.
   - `fitting/mod.rs` — `ScanResult<P>`, `FitTermination`. The scan returns the first
     and the alternation reports the second; both are path-independent, which is why
     `ScanResult` is generic over the noise parameters.
   - `generic/mod.rs` — `SampleRates`, `GenericSampleParameters`, `CoupledFit`, and the
     four floors. All SNP/indel.
   - `generic/runs.rs` — `RunsModelStarts`, `RunsModelFit`, `StartOutcome`.
2. **`homozygous_non_reference_rate()` returns `Option`**, where the architecture writes
   a bare `GenotypeFrequency`. `by_alt_copies` is a public vector, so an empty one is
   representable and the bare return would have to panic. The accessor answers `None`
   for the same reason its sibling does — "this set does not have one" — rather than
   failing.
3. **`Ploidy` gained a `Display` impl** in `types.rs`. `GenotypeFrequenciesNotFittable`
   names a ploidy in its message, so the type has to render; it renders the bare number,
   because the message supplies the word.
4. ~~**`SampleRates`' simplex invariant is documented, not enforced.**~~ **Withdrawn
   after review, and it was half wrong on the half that mattered.** Deferring the *sum*
   check was defensible — no accessor reads the sum, so an unchecked sum could not
   produce a wrong number until Milestone E. Deferring the *length* check was not:
   `homozygous_non_reference_rate()` returns the **last** entry, so a ploidy-2 set
   holding one entry handed back the homozygous-*reference* rate — near 1.0 — under the
   homozygous-*non-reference* name, where the truth is near 0.001. What landed instead is
   a checked `try_new` with private fields, and the error-variant objection is answered by
   the `.expect()` pattern the four constrained scalars already use.

## 3. Changes made

**`parameter_estimation/mod.rs`** — `Provenance` (four variants, and the doc says the
failure the type exists to prevent is a consumer that treats them alike);
`Estimate<T>`; `ParameterEstimationError` with the four variants and a `#[from]`
`DomainError`.

**`fitting/mod.rs`** — `ScanResult<P>`, generic over the noise parameters so the STR
path's three stutter parameters fit the same type, carrying `argmax_at_ladder_end`; and
`FitTermination`.

**`generic/mod.rs`** — the four floors, `SampleRates` with `observed_heterozygosity()`
and `homozygous_non_reference_rate()`, `GenericSampleParameters`, `CoupledFit`.

**`generic/runs.rs`** — `RunsModelStarts` with its `Default` (three separations ×
three inside fractions), `StartOutcome`, `RunsModelFit`.

**`types.rs`** — `impl Display for Ploidy`.

## 4. Tests added

Fourteen, across four files. The ones with teeth:

- `each_fitting_failure_names_the_sample_and_the_number_that_was_too_small` — the
  plan's stated A6 oracle. Each message must carry the sample, the number that fell
  short, and the floor it fell short of; these are read out of a log on a cohort run of
  hundreds of samples.
- `the_unseparated_states_message_says_it_is_not_an_inbreeding_coefficient_of_zero` —
  the one message that has to say what it is *not*. An outcrosser and a failed search
  leave identical fitted values, so a reader who takes this for `F = 0` has been handed
  a confident wrong number, which is the whole reason the variant exists.
- `heterozygosity_is_absent_above_diploidy_where_the_homozygous_rate_is_not` — the
  plan's stated A5 oracle, plus the frequencies summing to one.
- `a_haploid_region_has_two_classes_and_no_heterozygosity` — the *other* side of the
  boundary. **Corrected after review:** this first said a `>=` would let a haploid
  answer. It would not — `1 >= 2` is false. Measured, `>= 2` lets a *tetraploid* answer
  and the tetraploid test catches it; `<= 2` lets a haploid answer, and only this test
  does. The test is right; the reason given for it was not.
- `the_default_starts_disagree_about_the_state_separation_not_only_about_f` — the
  property `RunsModelStarts` exists for. Starts differing only in the assumed inside
  fraction are not a spread: they miss a genome whose states sit close together in the
  same way, and "keep the best-scoring" then has nothing better to pick. The measured
  cost is `F` = 0.2634 against a converged, silent `F` = 0.0000 on the same genome.
- `the_default_error_rate_is_a_constructible_rate` — the defaulted rung of the error
  rate's fallback ladder cannot be taken at all if the constant is not a value the type
  accepts.

## 5. Validation results

All in the container.

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo doc --no-deps --lib` | 12 unresolved links, all pre-existing; none in this step |
| `cargo test --lib ng::parameter_estimation` | 34 passed, 0 failed |
| `cargo test --all-targets --all-features` | 2,934 passed, 1 failed (pre-existing), 5 ignored |

The one failure remains `ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes`, confirmed pre-existing at `HEAD`.

## 5a. Review outcome

Three agents, eight categories, each in its own worktree: 0 Blocker, **6 Major**, 10
Minor, 4 Nit; 18 applied, 4 deferred, 0 disputed. See the
[review](../reviews/ng_parameter_prepass_generic_a5a6_2026-08-06.md) and the
[fixes applied](../reviews/fixes_applied_2026-08-06_v3.md).

**Three of the six Majors are wrong numbers rather than wrong code**, which is where the
defects live on a milestone whose deliverable is declarations. One is the `SampleRates`
inversion above; one is "29% covered by runs" against a realised 0.2629 — copied from
`spec` §6.5 and `arch` §5.3, **which carry the same slip and are not corrected here**;
and one is in this report, above.

Eleven other numeric claims in these files were checked against the research note, both
specs and the architecture, and all eleven held exactly.

After fixes the module holds 40 tests (was 34) and the suite is 2,941 passed / 1 failed
(pre-existing) / 5 ignored.

## 6. Tradeoffs and follow-ups

- **`ParameterEstimationError` lives in `mod.rs` while the floors it quotes live in
  `generic/`.** The error is the step's — `#[non_exhaustive]`, and the STR path adds
  variants — while the floors are the SNP/indel path's measured numbers. The parent
  imports two constants from the child for its messages, which is the direction that
  costs nothing.
- **The messages interpolate the floors by name**, so a change to `MIN_SITES_TO_FIT`
  changes the message with no second edit.
- **Nothing constructs a `GenericSampleParameters` yet.** Milestone E fills it; F1 is
  the entry point that returns it.
