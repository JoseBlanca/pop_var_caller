# Each library's own stutter, step H1 — the fit keeps each library's reads apart

**Date:** 2026-10-10. **Plan:** [str_library_stutter.md](../../implementation_plans/str_library_stutter.md)
step H1. **Spec:** [str_library_stutter.md](../../ng/spec/str_library_stutter.md) §2.1, §3.1 item 1,
§3.5. **Branch:** `stutter-multiplier`. The code review and the fixes it led to are §5 of this
report; the reviewers' own notes are not kept.

## 1. What was built

The repeat-tract fit used to add up, at every tract, the reads of all of a sample's read groups that
share a slippage group, so it could not score one library's reads differently from another's. **It
now keeps each library's reads apart and scores each library at its own level: the class rate times
the library's stutter multiplier.** Every multiplier is one in this step, so no number moves. The
multipliers are fitted from step H3 on.

- **`Library`** (`ssr_fit.rs`): one read group of the cohort, the slippage group whose class rates it
  reads, and its `stutter_multiplier`, which is one until the fit sets it.
- **`StratumEvidence::libraries`**: the cohort's libraries, in the cohort's read-group order. **One
  list a run, shared by every stratum** (`Arc<[Library]>`), so a later step that updates the
  multipliers hands every stratum the same new list and none can be left at an old multiplier.
  `SampleTractReads::by_library` (was `by_group`) holds a sample's reads per library, named by a
  `LibraryIndex` into that list — a type of its own, so it cannot be confused with the sample or the
  slippage group.
- **`library_level(class_level, multiplier)`** gives the level a library's reads are scored at:
  `min(m · L, max(L, MULTIPLIED_LEVEL_CEILING))`, with the ceiling at 0.95; a multiplier of one
  returns the class rate itself.
- **`gather_strata`** makes every read group of the cohort a library. A read group that
  `slippage_group_of` does not name is put in slippage group 0, which is what happened to its reads
  before; the log now counts such read groups.
- **Scoring (`ReadDistributionsByLevel`)** builds one read-length distribution for each slippage
  group and multiplier that some library is scored at, not one for each library. A sample with one
  library runs the loop the fit ran before. A sample with several libraries at the same level has
  their counts added together before the logarithm, as the old evidence did. With every multiplier
  one, every sum is the same as before, to the bit.
- **`StratumEvidence::check_libraries`**, called where a fit first takes a stratum, refuses by name a
  library in a slippage group the evidence does not declare, and a sample's reads naming a library
  past the list or out of order.
- **The subset of samples a large cohort's stratum is fitted on** still tops up each *slippage
  group* to eight samples with reads, not each library (spec §3.5).
- **The period pools** refuse strata whose library lists name different read groups or groups;
  different multipliers alone are not refused.
- **`draw_stratum`** (the drawn-stratum fixture) and the two examples that build strata by hand now
  make each sample a library of its own, in group 0 (`one_library_a_sample`).

## 2. Departures from the spec, both small

1. **The ceiling limits how far a multiplier can raise a level. It never lowers a class rate.** The
   spec writes the level as `min(m · L, 0.95)`. If the class rate `L` is above 0.95 — a curve may
   give up to 0.999 — that formula gives 0.95 even at a multiplier of one, so the fit would no longer
   match today's to the bit (spec §3.1 item 1) and calling would change for any such class. The
   implemented formula is the same as the spec's whenever `L` ≤ 0.95. Above it, it keeps the class
   rate, and a multiplier above one cannot lift the level further. **The spec's §2.1 wants amending
   to match; that waits for the owner.** Its §3.2 bound on the multiplier, `LEVEL_CEILING / maxₛ Lₛ`,
   is to be re-checked against this cap at step H2.
2. **The constant is named `MULTIPLIED_LEVEL_CEILING`, not `LEVEL_CEILING`.** The curves already
   have a `LEVEL_CEILING` (0.999, `slippage_curve.rs`).

## 3. How it was checked

- **`libraries_kept_apart_at_a_multiplier_of_one_fit_as_their_sum_to_the_bit`.** Six drawn strata
  of 8 samples are fitted twice: with each sample's reads split between two libraries (half of each
  bucket, rounded down, to one and the rest to the other), and with each sample's reads in one
  library. Every fitted number, standard error and walk record agrees to the bit, with all samples
  in one slippage group and with samples 5 to 7 in a second.
- **Scoring**: `a_read_distribution_is_one_slippage_group_at_one_multiplier` (which libraries share a
  distribution, and in what order); `tract_likelihoods_do_not_depend_on_the_order_a_samples_libraries_are_listed_in`;
  `a_samples_libraries_at_two_multipliers_are_two_rows` (two libraries of one group at multipliers 1
  and 2 score as two rows, not as their sum).
- **Gathering**: `gathering_keeps_each_read_group_a_library_in_its_declared_group` (`census_fit.rs`),
  on the fixture cohort of three read groups over two samples, gathered with each read group in its
  own slippage group.
- **Refusals**: a library past the declared groups, reads naming a library past the list, and a
  period pool of strata naming different libraries.
- **`a_library_level_is_the_class_rate_times_the_multiplier_below_the_ceiling`**, including a level
  that is not a number at a multiplier of one.
- **Before the review's fixes**: the fit's unit tests in release (75 passed, 4 ignored), the full
  debug suite (all green but the three known example failures, `ng_generic_loci_dump` × 2 and
  `ng_ssr_loci_dump` × 1), and the oracle (`scripts/promote_ng_oracle.sh`, four tomato accessions,
  one read group each, 7,451 records): every checksum matches the baseline. **After the fixes**, by
  the owner's lighter rule for a step that moves no number (2026-10-10): the tests above and the
  size test, in release, and clippy. The oracle and the full suite run again at the milestone.

## 4. Size and speed

- **Evidence size.** One library list per run, 16 bytes a read group, shared by every stratum and
  not counted in a stratum's `heap_bytes`.
- **Scoring speed.** Measured on one drawn stratum of 400 tracts, 200 samples at three reads, 13
  allele classes (the distributions built once, then every tract's likelihood table), best of seven,
  natively on the owner's Mac, this branch and the commit before it run alternately four times:
  0.096, 0.092, 0.080, 0.082 s against 0.096, 0.098, 0.089, 0.093 s, every run's checksum the same.
  The machine carried a load of about 12, so the figures are noisy; no slowdown is visible.

## 5. Review and fixes

Two reviewers, one on correctness (eleven deliberate defects planted; three went undetected by the
tests; each is now covered by a test built on the reviewer's probe, which failed under it) and one on design and speed.

| finding | outcome |
|---|---|
| scoring 1.75 times slower on one-library samples (0.31 → 0.55 s, reviewer's measurement) | fixed: one-library samples run the old loop; §4 |
| a sample's rows in different groups, order untested | fixed: test added |
| `gather_strata`'s library list and read keying untested | fixed: test added in `census_fit.rs` |
| library list copied into every stratum | fixed: one shared list |
| library index a bare `u32` | fixed: `LibraryIndex` |
| `groups` and `libraries` unchecked against each other | fixed: `check_libraries` at the fit's entry |
| a level that is not a number turned into 0.95 | fixed: multiplier one returns the level; a debug check on the multiplier |
| pool check compared multipliers as floats | fixed: compares read group and slippage group |
| name `ScoringTables` taken by another module; `levels` held keys | fixed: `ReadDistributionsByLevel`, `group_and_multiplier` |
| quadratic table lookup; scratch vectors reallocated | fixed: binary search; vectors reused |
| undeclared read groups silently in group 0 | fixed: counted in the log, stated in the docs |
| two libraries at two multipliers untested; pool refusal untested | fixed: tests added |
| the multiplier unchecked as a type | deferred to H2, where the first multiplier other than one is made |
| `library_level` belongs beside the stutter model calling reads | deferred to H5, where calling first uses it |
| spec §2.1's formula differs from the code | for the owner (§2) |
