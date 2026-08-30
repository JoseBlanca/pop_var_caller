# Fixes applied — ng parameters file, C2 first half

**Date:** 2026-08-30
**Review:** [ng_parameters_file_c2a_2026-08-30.md](ng_parameters_file_c2a_2026-08-30.md)
**Outcome:** all four Majors and twelve of fourteen Minors applied; two Minors recorded with a
reason.

---

## Applied

**Correctness**

- Four read-group-keyed tables and two per-sample tables now have to cover the axis they are keyed
  by, with no gap and no duplicate — the hole that let a deleted calibration row become a silently
  defaulted one.
- *Not measured* is a disjunction, matching `ContaminationView::was_measured`; the row with zero
  markers and non-zero reads is now refused, and it is the worse of the two.
- The contamination share is half-open at one, matching the consumer's own assert rather than
  reaching it as a panic.
- Eleven floats in the repeat-tract provenance are checked for being numbers.

**Messages** — every path is now a literal substring of the file; the read-group gap names the
missing id rather than a healthy row; three messages say what to do rather than stating a fact.

**Prose** — four false claims corrected, two of them mechanisms: why no line number is given, and
who calls `validate`.

**Tests** — 21 became 24, and the path test now checks every segment rather than the last, which is
what would have caught the wrong-key defect.

## Recorded, not applied

- **The saturation marker is checked on counts and not on the nine other integers the writer can
  saturate.** Those are provenance — repeat counts, cell counts, stratum counts — rather than
  evidence behind a number, and the header now says which are checked. A saturated
  `reference_repeats` is a lost number and would be worth catching; it is not worth another nine
  call sites in this step.
- **A warranted value may carry `observations = 0`**, where the shape's doc says a count is absent
  rather than zero. Its consequence is a wrong number in a run's report, not a wrong score, and it
  is the only member of that family C2 does not refuse.

## Validation after the fixes

- `cargo test --lib ng::calling::parameters_file` — **103 passed, 0 failed, 2 ignored**.
- `cargo test --all-features --lib --tests` — **5,026 lib tests plus every integration binary, 0
  failed**.
- clippy and `fmt --check` clean.
- Mutation sweep re-run: nine mutations, nine killed.
