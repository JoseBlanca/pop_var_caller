# Review — window coverage A3: a depth axis fitted to the sample

**Date:** 2026-09-06
**Reviewed:** commit `8b782af1` on `ng-window-coverage` (plan step A3 of
[window_coverage.md](../../ng/impl_plan/window_coverage.md))
**Impl report:** [ng_window_coverage_a3_2026-09-06.md](../implementations/ng_window_coverage_a3_2026-09-06.md)
**Per-category files (audit trail):** `tmp/review_2026-09-06_window-coverage-a3/`

## 1. Scope

The A3 diff — the depth bin width becoming a per-sample fit — across
`src/ng/window_coverage/{mod.rs,accumulator.rs,production_parity.rs}`. Two agents in isolated
worktrees split the categories: one on `reliability`, one on `naming`, `idiomatic`, `defaults`,
`errors`, `refactor_safety`, `smells` and `module_structure`.

Out of scope: `src/sample_summary/` (frozen; read as the oracle); the four checks red on `main`
before this branch.

## 2. Verdict

**Approve with changes.** Both agents independently found the same defect, and it is the exact
failure the plan singles out for this step — "its failure is silent — a width off by a factor is
a histogram whose mode sits in the wrong bin, and every fit anchors on it".

## 3. Top three

1. **A fit that failed part-way through a sample silently retried, and dropped the windows it
   failed on.** The fit took the held-back windows out of the accumulator, found no positive
   median, and returned — leaving the width unset *and* the buffer empty, so the next windows
   started a second scale sample. The sample then got a histogram whose width was fitted from a
   later stretch of its genome, with the earlier windows folded nowhere and the count agreeing
   with the cells. Measured on a six-window fixture with depths `0,0,0,80,30,30`: six windows
   emitted, none absent, **`windows_folded` = 2**. Two doc comments promised the opposite.
   **Fixed** by replacing the width's `Option<f64>` plus a separate buffer with a three-state
   type — collecting, fitted, or **unfittable and latched** — so a sample that cannot be scaled
   stays unscaled. The regression test is the fixture above.
2. **The depth axis was not protected against a contig change.** Adding a line to `reset_contig`
   that cleared the fitted width left all 42 tests green; a probe showed a width of 2.5 becoming
   10 when a deeper second contig followed — one histogram whose rows are cut on two different
   axes. No fixture fitted a width before a contig boundary. **Fixed** with
   `the_fitted_depth_width_survives_a_contig_change`, and the field's doc now says why the width
   is per sample and never reset.
3. **The differential's narrowing threw out two properties that had nothing to do with the depth
   axis** — that a stream of real windows produces a histogram at all, and that every window
   emitted reaches it. Above 25 windows nothing guarded either: a mutation returning no histogram
   for any sample with more than 100 scale windows, and one that stopped folding after the 100th
   window, both survived the whole suite. **Fixed**: the differential asserts both again, which
   costs nothing and does not depend on how either side cuts its depth axis.

## 4. What's good

- **The transcription still holds.** The differential compares 447,581 window means and GC
  fractions bit for bit with production, the same count as before A3 — so the fit changed nothing
  about what a window reports.
- **Refactor safety is what surfaced A3's new fields at all.** The configuration has no
  `Default`, so each new field is a compile error at every literal, and `finish` destructures it
  exhaustively with named ignores rather than `..`.
- **The histogram's stated size checks out exactly**: 50 × 401 × 4 bytes is 80.2 kB.

## 5. Findings

### Blocker

**B1 — a failed fit retried and lost its windows.** *Categories: errors, reliability
(convergent).* Described above. **Applied.**

### Major

**M1 — the fitted width was unguarded against a contig reset.** *Category: reliability.*
Described above. **Applied.**

**M2 — the differential dropped two axis-independent properties.** *Category: reliability.*
Described above. **Applied.**

**M3 — the window's two means travelled as a bare `(f64, f64)`, read by `.1`.** *Category:
naming.* Two same-primitive domain scalars, four positional hand-offs, and the one that decides
the sample's whole depth axis identified its field by index. Transposing them is silent in every
sense the module has: the emitted pair is built separately, so the per-locus half stays right and
the differential stays green, while the width is fitted from a median **GC fraction** and comes
out hundreds of times too small. **Applied:** a named `WindowMeans` pair, and
`median_depth_returns_the_upper_middle_on_an_even_count`, which asserts that the GC fraction is
never the key even when it orders the pairs the other way.

### Minor

**Mi1 — `depth_range_in_medians` was typed `u32`, a count for a multiple.** *Category:
idiomatic.* Its only use converts it straight to floating point, and `>= 1` additionally forbade
a range below one median — a legitimate setting for a sample whose median sits above its mode.
Plan step D2 is what settles this number, and `u32` said in advance that D2's answer had to be a
whole number. **Applied:** `f64`, with `is_finite() && > 0`. This also re-armed a warning A3 had
left guarding a field it had just deleted — `assert_valid`'s doc still described a `NaN`
`depth_bin_width` after `depth_bin_width` left the configuration, so at `8b782af1` no field could
be `NaN` at all. Now the `NaN` case is real again (a `NaN` multiple makes the fitted width `NaN`,
which the fit reads as "this sample cannot be scaled", so a misconfigured run would cost every
sample its histogram silently) and `new_panics_on_a_nan_depth_range` pins it.

**Mi2 — `DEPTH_RANGE_IN_MEDIANS` did not say it was provisional** where it is defined, while the
other three soft constants did. It is the constant most likely to move: D2 measures the overflow
fraction against the fit's rejection guard, and that fraction is a direct function of this number
and nothing else. **Applied**, along with the configuration doc that counted ng's additions as
two when there are three, and the exhaustive-destructure comment that said "neither" of three.

**Mi3 — one quantity, two vocabularies.** The fitted value was `depth_bin_width` where it is
stored and read, and "the scale" where it is fitted (`windows_awaiting_the_scale`,
`fit_the_depth_scale_and_fold_what_was_held_back`) — one letter from `depth_scale_windows`, which
is a count of windows and not a scale at all. **Applied:** one vocabulary, and the method is now
`fit_depth_bin_width_and_fold_held_back`.

**Mi4 — the scale sample's boundary was untested.** `>=` slipping to `>` changed a fitted width
from 2.5 to 1.0 with the suite green. **Applied:**
`the_scale_sample_is_exactly_as_many_windows_as_it_says`.

**Mi5 — `median_depth` had no direct test**, and the held-back fold's only test held depth
uniform, so a defect present in both of its arms cancelled. **Applied:** a direct test over
empty, one, even, odd, all-equal and transposed inputs, and
`a_held_back_window_lands_in_the_same_cell_as_one_folded_immediately`, whose depths differ and
whose cells are named by hand.

**Mi6 — `sliding_uniform_all_one_cell` asserted a shape, not a cell.** Twenty identical windows
land in one cell under any positive width and any in-range index formula. Before A3 the width was
configured and the cell was computable; the fit made it a function of the data and the test was
not updated. **Applied:** the width and the cell index are both named.

**Mi7 — the memory note was right per window and wrong in total.** Sixteen bytes each is correct,
but the list grows by doubling, so at 10,000 windows its capacity reaches 16,384 and the
transient peaks at **262 kB, not 160**. **Applied** to the code's doc; the spec quotes 120 kB and
is owed a correction — see §7.

**Mi8 — `finish`'s `Option` reports three different silences as one.** *Category: errors,
question 4.* No window finalised, every window under the floor, and a non-positive median are
three different bug reports, and at the project's declared hardest corner — one low-coverage
sample — the middle one is the likely case and is a reading on `MIN_WINDOW_POSITIONS`, which plan
step D1 exists to settle. The reviewer implemented a four-variant `SampleHistogram`.
**Deferred and raised with the owner at Checkpoint A**: spec §3.6 fixes `finish`'s signature and
§3.5 its prose, so this is a design change rather than an implementation choice.

### Nits

Applied: the stale A2 doc on `an_absent_window_equals_itself_and_signed_zeroes_differ`, which
still said nothing in the module produced an absent pair; the module doc that called the depth
width "still to come" after A3 landed it; and a sentence on the early/late fold test naming what
a uniform depth cannot see.

Considered and not taken: renaming `median_depth`, a noun-named free function that reorders its
argument (both are documented at its definition); splitting `WindowCoverageConfig`, which at six
fields is still under the smell's threshold and whose fields are all chosen at one site and
validated together; and pulling the cell scheme into a `DepthAxis { bins, width }` built at the
fit, which the reviewer recommends only if a seventh field arrives.

## 6. Verification

- **Mutation testing.** The reliability agent ran 16 mutations against A3 as committed: 11
  caught, 5 survived, 2 of those changing no behaviour. The design agent ran its own and found
  the same Blocker independently.
- **Five mutations re-run by the orchestrator on the fixed tree**, each caught by the test
  written for it: the unfittable latch removed (fails `a_failed_fit_does_not_start_a_second_scale_sample`
  alone); the width cleared at a contig change (fails `the_fitted_depth_width_survives_a_contig_change`,
  and the differential's restored assertion); the scale sample's `>=` slipped to `>` (fails
  `the_scale_sample_is_exactly_as_many_windows_as_it_says` alone); the median taken on the GC
  fraction (fails `median_depth_returns_the_upper_middle_on_an_even_count` and 12 others); and
  the held-back windows folded against the median rather than the width (fails
  `a_held_back_window_lands_in_the_same_cell_as_one_folded_immediately` and 4 others).
- `cargo test --lib --all-features` — **6,323 passed, 0 failed, 15 ignored**; 48 in this module.
- `cargo clippy --lib --all-features` — no diagnostic in `src/ng/window_coverage/`.
  `rustfmt --check` clean.
- **The calling oracle, re-run after the fixes** — six tomato accessions called from the CRAMs and
  from stored files: 2,311 records, sha256 `84ad19c2…` on both sides, identical to the baseline
  taken before A1.

## 7. For the owner, at Checkpoint A

- **`finish`'s return shape** (Mi8) — recommended, and it changes spec §3.5 and §3.6.
- **Spec §3.4's 120 kB transient** is 262 kB. Its budget paragraph still closes — 123 + 80 + 262
  = 465 kB against 500 — but with about 35 kB of headroom rather than the ~180 the quoted figure
  implies. Plan step D3 measures it; worth correcting before it is quoted again.
- Carried from the earlier steps: spec §3.4's "Every finalised window is folded" contradicts the
  code; a `windows_under_the_floor` counter on the histogram; and spec §5's 12-byte deque entry,
  which is 24 bytes.
