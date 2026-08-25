# Fix Application Report: ng_calling_loop_c1_2026-08-25.md

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`
**Review:** [ng_calling_loop_c1_2026-08-25.md](ng_calling_loop_c1_2026-08-25.md)

## 1. Executive summary

**0 Blockers, 8 Majors, 5 Minors. All 8 Majors fixed** — five by changing code or tests, three by
correcting claims. Tests went from three to **seven**; the library target went 4,599 → **4,603**.

### The one behaviour change

**Two release-held shape checks were lifted out of the `match`**, so they hold on both arms.
`prior_row`'s width was `PriorRow::new`'s to check and `sample_expected_copies`' was
`fill_sample_concentration`'s, and a flat pass enters neither — so a mis-shaped buffer made the
seeded arm panic in release and made the flat arm return a wrong posterior **in silence**. Two
`#[should_panic]` tests now reach them on the flat arm.

### Two claims corrected, both mine

1. **The window was reported as though cohort size were not an axis.** The sweep's own output
   shows the two starts diverging at 0.5 nats for 20 and 63 samples; the comment said there was
   "nothing to lose" there and that the effect bit near 1 nat "and nowhere else". The comment now
   gives both axes and says plainly that the earlier version was the failure `CLAUDE.md` names.
2. **The rejected design's failure mode was wrong in both halves.** It panics on this function's
   own release-held own-copies check, **in release as well as debug** — not on the cohort check,
   and not silently. The bare-seed fall-through is real only when the cohort row alone is `NaN`.
   `PassPrior`'s doc now separates the two cases and gives the probe for each.

### The claim that was corrected by measurement rather than by argument

**The trap is a delay, not a different answer, and the test now says so.** Measured at 63
samples: the seeded start sits at 0.151 alternative copies per sample at six passes, **0.633 at
nine, where it flips to heterozygous**, and both starts reach 0.767332 by thirty. At three
samples the flip is between ten and sixteen. A rare-variant fixture — a handful of carriers among
60 firmly homozygous-reference samples, swept over carrier count and advantage at 50 passes —
showed the two starts agreeing in **every** cell.

The test is renamed
`a_flat_first_pass_finds_the_variant_at_once_where_a_seeded_one_takes_nine_passes`, carries the
pass-by-pass table, and **asserts the agreement at thirty passes** so the "delay, not divergence"
fact is pinned rather than stated. `spec/calling_em_loop.md` §3's stronger claim is flagged for
the owner in the review's §6.

### Validation

| command | result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `4603 passed; 0 failed; 14 ignored` |
| `cargo test --release --lib ng::calling --all-features` | `557 passed; 0 failed` |

## 2. Per-finding log

### M1 — the flat arm bypassed two shape checks — **Fixed**

Both lifted above the `match`, with a comment recording the measured silent-wrong-answer they
prevent. Tests: `a_prior_row_of_the_wrong_width_is_refused_on_a_flat_pass` and
`an_expected_copies_row_of_the_wrong_width_is_refused_on_a_flat_pass`, both reached with
hand-built bundles because the scratch cannot produce a short row.

### M2 — nothing tested in release that the flat pass leaves the cohort alone — **Fixed**

`the_flat_pass_touches_neither_the_cohort_summary_nor_the_concentration` asserts *positively* —
after a flat pass both the cohort row and the concentration still hold the `NaN` sentinel — so it
holds in both profiles, where the previous cover was a `debug_assert` one module away.

### M3 — the trap test's stated claim was refuted by its own fixture — **Fixed**

Above.

### M4 — a silent sample votes for a 50% allele frequency — **Pinned, and raised**

`a_sample_with_no_reads_contributes_a_full_alternative_copy_to_the_flat_pass` records the
behaviour and the number (`1.0` copies of each allele). **Not changed**, because neither §3 nor
§7 says what should happen and the choice would shape the design: §7's *"the prior decides it
alone"* has no prior to appeal to on a flat pass. Raised for the owner.

### M5–M8 and the Minors — **Fixed**

- `reads_the_cohort`'s doc claimed to be "the one place the two arms are told apart"; there are
  two, and it now says so.
- **`CohortFixture` had been inserted between `run_passes`'s doc comment and `run_passes`**, so
  29 lines of doc documented the wrong item and `run_passes` had none. Moved back, and
  `CohortFixture` given a doc of its own.
- The window and failure-mode corrections above.

## 3. Not applied

- **Renaming `PassPrior` to something per-sample.** It is a per-*sample* value under a per-*pass*
  name, which a reviewer demonstrated by writing D1's loops — the per-sample form compiles and
  the hoist the name invites fails with `E0425`. The name reads correctly at the call site and in
  the spec's own vocabulary (*"the first pass runs on the reads alone"*), and D1 is where a loop
  will show whether the invitation is a real hazard. Recorded rather than pre-empted.
- **`Debug` on `PassPrior`.** `GenotypePriorModel` has no `Debug` supertrait — the reason
  `name()` exists at all — so deriving it does not compile without one.

## 4. Carried forward

1. **Spec §3's "converges to no-variant"** — not reproduced; the owner's to settle.
2. **A silent sample's flat-pass vote** — a design gap between §3 and §7.
3. **The flat start is a property of three samples and up on this fixture**; at two samples it
   does not hold the heterozygote either, and at one spec §7's claim holds bitwise.
