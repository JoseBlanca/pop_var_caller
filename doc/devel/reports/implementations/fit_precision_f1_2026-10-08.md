# Fit precision, step F1 — own-fit errors beside blends, and every stratum counted in the log

**Date:** 2026-10-08 · **Branch:** `fit-precision-followup` · **Plan:**
[fit_precision.md](../../implementation_plans/fit_precision.md), "After the kimura run — decided" ·
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §5.2, amended by the owner on 2026-10-08.

## What changed, and why

The owner's kimura run (2,169 samples) wrote no own-fit standard error for 31 of the 34 repeat-tract strata
fitted on their own tracts, because their slippage level is a blend of the stratum's own fit with its
period's curve, and spec §5.2 withheld an error from any blended number. In the 22 one-base strata among
them the curve's weight in the blend was 2 × 10⁻⁵ to 1 × 10⁻⁴: the numbers were, in effect, the own fit. The
same run's log accounted for 61 of its 141 strata; the other 80 were refused and no line said so.

- **A slippage number keeps the own fit's error wherever the own fit went into it** — the own fit whole, or a
  blend of it with a curve. Only a number taken whole from the curve has none. The error is still the own
  fit's, from before the blend; the blend's own is not computed. The rule is one exhaustive match a source:
  `LevelSource::holds_the_own_fit` and `ShareSource::holds_the_own_fit`, used by
  `OwnFitStandardErrors::of_the_numbers_emitted` (`stratum_fits.rs`).
- **The file's reader accepts an error beside a blend** and refuses one only beside a number taken whole
  from the curve (`validate.rs`, `a_level_the_own_fit_went_into` and `a_share_the_own_fit_went_into`, both
  exhaustive).
- **The file's note on these errors** says the error beside a blend is the own fit's, and that the blend's
  `curve_weight` says how much of the number is the curve's: near 0 the error describes the number written,
  near 1 mostly the own fit the curve outweighed.
- **One new log line** once every stratum has its answer (`strata_outcomes_summary`, `ssr_fit.rs`):
  `strata: N in all; F fitted on their own tracts, D furnished from their period's curves, R refused — A
  with no read spanning a tract, B with fewer tracts with reads than the K a fit needs, holding x to y`.

## What moved, measured

**On four tomato accessions** (`scripts/promote_ng_oracle.sh`, 20 regions, the baseline's slice, 7,451
records): only the two parameters files' checksums move. Every VCF checksum, including the calls made with the
fit, is unchanged.

- The fitted file gains 36 own-fit errors — 14 levels, 12 shorter shares, 10 fall-offs, every one beside a
  blend. Its 44 blended numbers have curve weights from 0.05 to 0.9999, median 0.6: at four samples and about
  three reads a position the curve is often most of the number, the opposite of kimura. 8 blended numbers
  stay without an error because their own fit gave none.
- **Proof that nothing else moved** (`tmp/f1_explain_oracle.py`): with the old note put back and those 36 keys
  deleted, the fitted file hashes to the recorded `3338658c…`; the defaults run's file needs only the note put
  back to hash to `b7c0a262…`.
- The new log line on that cohort: `strata: 68 in all; 15 fitted on their own tracts, 22 furnished from their
  period's curves, 31 refused — 19 with no read spanning a tract, 12 with fewer tracts with reads than the 8 a
  fit needs, holding 1 to 5`.

**The three cross-platform fixtures** write no slippage rows, so they gain no key; each file's checksum moves
by the note alone. Putting the old note back reproduces each recorded checksum exactly
(`tmp/f1_explain_digests.py`). The calls' checksums do not move. Re-recorded in the Linux container and
checked on macOS (arm64).

## Tests

- `own_fit_errors_are_kept_for_every_number_the_own_fit_went_into` (`stratum_fits.rs`): three slippage
  groups of one fitted stratum, so each of the three numbers is seen from each of the three sources — the own
  fit, a blend, the curve whole. Putting the old rule back for the fall-off alone fails it (measured: the
  blended fall-off's error comes back `None`).
- `every_key_version_two_added_is_refused_where_it_means_nothing` (`validate.rs`): an error beside a level or
  a shorter share taken whole from the curve is refused; a blended fall-off with an error is accepted; the
  fixture itself carries errors beside a blended level and a blended shorter share, and passes.
- `the_strata_outcomes_summary_accounts_for_every_stratum` (`ssr_fit.rs`): the whole line, for seven strata of
  every kind, for one stratum below the floor, and for none refused.
- The golden files `every_shape.toml` and `every_shape_as_written.toml` gain the two blend errors and the new
  note.

## Deviations from the plan

- The review fanned out to two reviewers covering all nine categories between them, not one per category: the
  diff is about 150 lines of code and each reviewer's worktree costs a full build.
- The spec's table row was edited in place, with the old wording quoted in the amendment beside it, as the
  spec does at its other amendments; the plan's step text asked only for the amendment to be recorded.
