# Fixes applied — ng parameters file, A2

**Date:** 2026-08-28
**Review:** [ng_parameters_file_a2_2026-08-28.md](ng_parameters_file_a2_2026-08-28.md)
**Code:** [src/ng/calling/parameters_file/](../../../../src/ng/calling/parameters_file/)

---

## 1. What changed

Every Major and Minor was applied. The shape gained a fifth member and a second type:

- **`repeat_tract_outlier_weight` moved into `WarrantedValue`** (M1), so spec §8's three honest
  defaults are all marked. Five numbers now carry the shape where four did.
- **`EvidenceCount` replaced the bare `u64` count** (M2, m2), so the unit reaches the file:
  `observations = { reads = 812344 }`, `{ covered_positions = 180600412 }`,
  `{ bases_compared = 40122 }`.
- **The two wrong unit claims were corrected** (M2): an inbreeding coefficient is fitted over
  covered reference positions, not windows; a repeat-tract substitution rate over bases compared,
  not reads.
- **The wrong seam citation was corrected** (m1): `RunParameters::assemble` is what reduces the
  pre-pass's `Estimate<InbreedingF>` to a bare `Vec<InbreedingF>`.
- **The calibration's count is no longer asserted absent** (M3). Spec §3.3 asks for it, and it
  exists on the `Estimate<ErrorRate>` the same seam drops — so the fixture carries it and the
  field's doc says where B1 must get it.
- **The prose stopped calling five non-`Warrant` mechanisms "warrants"** (M4), and each keeps its
  own settled word: a contamination fraction *stands on* its evidence counts, the ordinary-site
  prior *carries a rung*, a slippage number *carries its origin*, a length spectrum is *placed
  rather than annotated*, and a read group's slippage group is *declared, not estimated*.
- **The two prose nits**: "the ordinary-site prior's two concentrations" for "the pair", and "a
  level interpolated off a curve" for "one read off a curve".

## 2. The one that is a fix forward into A1, not an A2 defect

M3 traces to a claim A1 shipped. A1's report states that the base-quality calibration has no
observation count, on the evidence of `ReadGroupCalibration` — which is the calling-side view, a
multiplier and a provenance. **Spec §3.3 asks for the count by name**, and the fit produces it: it
is on the `Estimate<ErrorRate>` that `RunParameters::assemble` reads and does not store. So A1
handed A2 a requirement that was wrong, A2 built to it, and A2's test then pinned it. The A1
report's §3 is left as it was written — the record of what was believed at the time — and this
paragraph is the correction.

## 3. What the fixes cost, and one deviation recorded

**`EvidenceCount` was not the reviewer's proposal.** `naming` proposed the opposite: drop the count
from the shape entirely, leaving `WarrantedValue { value, warrant }` with no optional key at all,
and put a unit-naming count on each row (`reads_behind_the_rate`, `windows_behind_the_coefficient`).
The argument for it was that the count is not in fact shared. **That argument weakened when
`reliability` established that the calibration does have a count**: three of the five warranted
numbers carry one, not two of four. Splitting a value from its count is the thing spec §2's spine
exists to prevent, so the count stays inside the shape and the unit rides on it instead. Both
reviewers' real complaint — that the unit was invisible in the artefact — is answered either way.

**The cost is a nesting level**, and under `serde`'s array-of-tables rendering that reads badly:
four header lines for one calibration row. Under step B2's inline writer it is one line a row.

## 4. Findings table

| # | severity | status | note |
|---|---|---|---|
| M1 | Major | **Applied** | `repeat_tract_outlier_weight: WarrantedValue`, with its two reachable states in the field's doc |
| M2 | Major | **Applied with adaptation** | units corrected, and moved out of rustdoc into the file as `EvidenceCount` |
| M3 | Major | **Applied** | the assertion inverted; the fixture carries the count; the field's doc names B1's source |
| M4 | Major | **Applied** | the reviewer's re-wording taken almost verbatim, extended to a fifth quantity |
| m1 | Minor | **Applied** | the seam is `RunParameters::assemble` |
| m2 | Minor | **Applied with adaptation** | see §3 — the unit is named, but on the count rather than on the row |
| m3 | Minor | **Applied** as documentation | `serde` gives no header to a section whose fields are all tables; both `[repeat_tracts]` and `[stated_constants]` have vanished. The module doc now says to read the golden file as a record of key names, not as the artefact |
| m4 | Minor | **Applied**, all three | the path list replaces the row count; the doc says which assertion catches which case; the stated-concentration test re-points at what a reader recovers |
| Nits | Nit | **Applied**, both | |

## 5. Validation

All in the container.

| command | result |
|---|---|
| `cargo fmt -- --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | **12 passed, 0 failed, 1 ignored** |
| `cargo test --lib` | **4,932 passed, 0 failed, 12 ignored** |
| `cargo doc --no-deps` | zero unresolved links in this module |
