# ng step 4, the STR path — E3: borrowing, against two floors

*Implementation report, 2026-08-12. Step E3 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — two agents; the reliability agent ran 38 mutations of which 15
survived, and I re-ran three decisive ones after the fixes. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.3 and §4.5,
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §4.1.*

## What the step is

`resolve_slippage` turns the searched fits into what each stratum will actually report, against
two floors that a stratum can fail independently:

- **fewer than `MIN_LOCI_TO_FIT` loci, or nothing its model could score** — the whole model is
  borrowed from neighbouring repeat counts, `Provenance::Borrowed`, with `fitted_over` naming the
  lenders;
- **fewer than `MIN_SLIPPED_READS_TO_FIT_SHARES` reads that moved** — the level is kept and only
  the direction share and the fall-off are borrowed, named in `shares_fitted_over`.

A period with no fittable stratum at all raises `NoFittableStratumAtPeriod` rather than defaulting.

## Recorded implementation choices

1. **Between two lenders the level is interpolated in the logarithm and the shares linearly, both
   weighted by distance in repeat counts.** A level rises about 1.3-fold per repeat count and
   spans orders of magnitude, so the multiplicative middle is what interpolates it; at the
   midpoint this is the plain geometric mean. Weighting matters when the lenders are not
   equidistant: with fitted strata at 5 and 12 repeats and everything between them thin, an
   unweighted geometric mean hands all six the same 0.0387 — 2.6 times too high at 6 repeats and
   2.6 times too low at 11.
2. **The level and the two shares are borrowed from two different sets of lenders.** Anyone fitted
   may lend a level; only a stratum that also cleared the second floor may lend shares.
3. **The second floor raises nothing when no stratum in a period clears it.** Spec §4.5 expects
   that wherever a whole period sits at the bottom of the repeat range, so it is the common case;
   each stratum then keeps what it measured. The state a reader needs — *its own, and nobody had
   better* — is `shares_fitted_over` empty beside a `slipped_reads` under the floor, derived
   rather than stored so there is one place the fact lives.

## What the review changed, and the first one would have stopped every real run

**A thin stratum stopped the whole sample before the borrowing ever ran.** The search fitted every
stratum and refused any whose four starts disagreed — and a two-locus stratum is exactly the one
whose starts cannot agree; measured, one locus gives a spread of 193 against a limit of 1.06. So a
genome with a single thin stratum would have died on a fit whose answer was about to be discarded.
**The locus floor now runs before the search**, in the walk, which also skips several hundred
searches on a real sample. The floor keeps one spelling, `thick_enough_to_fit`, used on both sides
of the seam because `resolve_slippage` is public and may be handed a map it did not build.

**A whole borrow laundered the shares the second floor had just rejected.** The lender for the
whole model was chosen on "was it fitted" alone, so a stratum whose own direction share of 0.99
was refused for standing on 40 moved reads would hand that same 0.99 to its thin neighbour as a
measurement. The two lender sets are now separate, and the two provenance lists can name different
strata.

**A borrowed number carried this stratum's evidence rather than its lenders'.** `observations` was
the borrower's own moved-read count — zero, for a stratum that had none — which reads as no
evidence at all. It is now the lenders' scored reads, summed, which is how the SNP/indel path's
own fallback settled the same question.

**The "no fittable stratum" message contradicted itself** when the blocker was that no read sat on
the whole-repeat grid: it said "no stratum with 1000 loci … its thickest holds 4000". It now
carries how many strata the period held and how many held a read that moved, and gives the two
remedies conditioned on those counts — widening the run helps one case and collects more of the
same in the other.

## The mutation run: fifteen survivors, one pattern

**Every fixture was alike where it mattered.** Slipped-read counts were 0, 40 or 9,000, so any
threshold between 41 and 8,999 behaved identically — and the locus floor's 1,000 sits inside that
window, so the second floor could read the wrong constant and pass. Every fixture used one motif
period and one library, so two of the three legs of the grouping filter cost nothing to delete.
No fixture ever handed a fit to a stratum under the locus floor, so the loop that clears one was
dead code. And `observations` and `loci` were read by no assertion at all.

**One survivor was an arithmetic coincidence in a test's own headline.** The stratum keeping its
level sat at 0.02 with neighbours at 0.01 and 0.04 — and `sqrt(0.01 × 0.04)` is 0.02 exactly in
binary floating point, so "kept, not borrowed" and "borrowed" produced bit-identical numbers. The
fixture now sits at 0.03.

Five new tests close them: a thin stratum arriving *with* a fit, the floor's boundary at 3,999
against 4,000, a period where nobody measured the shares, a borrow between distant lenders, and
borrowing across neither a period nor a library.

Re-run after the fixes rather than taken on the agent's word: deleting the thin-stratum floor,
pointing the shares floor at the locus constant, and dropping the library and period from the
grouping — all three now fail, at three different tests.

## Tests

Eleven for this step.

| test | what it pins |
|---|---|
| `a_thin_stratum_between_two_thick_ones_borrows_the_whole_model` | the two-sided borrow, the geometric level, the arithmetic shares, both lenders named, and the lenders' reads as the warrant |
| `a_thick_stratum_with_almost_no_moved_reads_keeps_its_level_and_borrows_its_shares` | the second floor, and the two provenance lists differing |
| `a_borrowed_model_does_not_take_shares_that_were_refused_next_door` | the level and the shares taken from different lenders |
| `a_borrow_between_distant_lenders_follows_the_gap_between_them` | 0.0147 at 6 repeats and 0.1019 at 11, between lenders at 5 and 12 |
| `a_thin_stratum_that_arrives_with_a_fit_borrows_anyway` | the locus floor applied to the map, not only inside the walk |
| `a_stratum_too_thin_to_fit_is_never_searched_and_so_never_refused` | the floor before the search — the run that would otherwise have stopped |
| `the_moved_read_floor_is_where_the_constant_says_and_the_boundary_keeps` | 3,999 borrows, 4,000 keeps |
| `a_period_where_nobody_measured_the_shares_keeps_them_and_says_nothing_was_borrowed` | the third state, and that it raises nothing |
| `a_period_with_no_fittable_stratum_errors_rather_than_defaulting` | the refusal, and both ways a period can fail |
| `a_stratum_at_the_end_of_the_range_borrows_from_the_side_it_has` | the one-sided borrow, taken whole |
| `borrowing_crosses_neither_a_period_nor_a_library` | and `borrowing_stays_inside_one_ploidy` for the third leg |

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` clean, and
`cargo test --lib --bins --tests --all-features` in the container: 3,458 → **3,470** lib tests, 0
failed, 11 ignored.
