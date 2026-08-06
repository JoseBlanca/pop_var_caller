# Fix Application Report: ng_parameter_prepass_generic_a5a6_2026-08-06.md

**Date:** 2026-08-06
**Source review:** `doc/devel/reports/reviews/ng_parameter_prepass_generic_a5a6_2026-08-06.md`
**Source state reviewed against:** `b9ef1e8`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 6 · Minors: 10 · Nits: 4

### Outcome totals
- Applied: 17 · Applied with adaptation: 1 · Deferred: 2 · Disputed: 0
- Already fixed / Failed validation / Blocked / Superseded / Awaiting answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib ng::parameter_estimation` → 0, **40 passed** (was 34)
- `cargo test --all-targets --all-features` → 101, **2,941 passed, 1 failed, 5 ignored** — the failure is the pre-existing `ng::locus_generation::pileup::parity` divergence
- `cargo doc --no-deps --lib` → 101, 12 unresolved links, all pre-existing
- Performance check → not applicable; nothing here is reachable from `benches/`

### Unresolved high-priority findings
None. Both deferrals are design-document changes outside this workflow's remit.

## 2. Findings table

| ID | Severity | Title | Final status |
|---|---|---|---|
| M1 | Major | `SampleRates` accessors read by dosage with no guarantee | Applied with adaptation |
| M2 | Major | Two documented measurements wrong | Applied (code); upstream deferred |
| M3 | Major | The impl report's stated reason for a test was wrong | Applied |
| M4 | Major | `Display for Ploidy` untested | Applied |
| M5 | Major | `RunsModelStarts::default()`'s fields swappable undetected | Applied |
| M6 | Major | `Domain` variant breaks three error rules | Deferred |
| Mi1 | Minor | `separations` is inverted | Applied |
| Mi2 | Minor | `Default` undocumented, values unstated | Applied |
| Mi3 | Minor | Floors lack citations and soft/fixed labels | Applied |
| Mi4 | Minor | SNP/indel floors hard-wired into a shared error | Applied |
| Mi5 | Minor | `MIN_WINDOWS_TO_FIT_INBREEDING` outside `runs.rs` | Applied |
| Mi6 | Minor | Raw `f64` where newtypes exist | Applied (`LogProb` only) |
| Mi7 | Minor | `starts_tried` ordering has no accessor | Applied |
| Mi8 | Minor | Floor assertions pass vacuously | Applied |
| Mi9 | Minor | Doc conflates the rate with a sum | Applied |
| Mi10 | Minor | `inside_het_floor` named a floor, documented a rate | Applied |
| N1 | Nit | `SmallVec` inline size contradicts its own doc | Applied |
| N2 | Nit | `#[non_exhaustive]` justified by an in-crate reason | Applied |
| N3 | Nit | `implied_f` reads as a formula symbol | Deferred |
| N4 | Nit | Ten remaining raw `f64`s in `runs.rs` | Deferred |

## 3. Questions asked and answers

None asked of the user. Two open questions from the review are recorded in `PROJECT_STATUS.md`.

## 4. Per-finding log

### M1 — `SampleRates` (Applied with adaptation)

- **Implementation.** `ploidy` and `by_alt_copies` are now private, behind
  `SampleRates::try_new(ploidy, by_alt_copies)`, which rejects a set that does not hold
  one frequency per dosage `0..=ploidy` or does not sum to one within a stated
  tolerance. `homozygous_non_reference_rate()` returns a bare `GenotypeFrequency` —
  there is no longer an empty case to guard. `observed_heterozygosity()` keeps its
  `Option`, which now has exactly one meaning: this genome does not have one.
  A `GenotypeFrequenciesOffSimplex { ploidy, entries, total }` variant carries the
  failure.
- **Adaptation, and it settles a disagreement between two reviewers.** `reliability`
  proposed a six-line `is_one_entry_per_dosage()` guard on both accessors, keeping the
  fields public — arguing that a checked constructor would need an error variant for a
  condition only our own arithmetic can produce, which was the commit's own objection.
  `idiomatic` proposed private fields and `try_new`. The constructor is what landed: the
  guard leaves both accessors answering `None` for a malformed set, and the doc defines
  `None` as "this genome does not have one of these" — so the guard would make the
  wrong-length case indistinguishable from the ploidy-4 case at every call site. The
  error-variant objection is answered by the same `.expect()` pattern the four
  constrained scalars already use: our own arithmetic being broken is exactly what a
  checked constructor plus `.expect()` is for.
- **Verification.** Removing the length check from `try_new` fails
  `a_set_with_the_wrong_number_of_entries_for_its_ploidy_is_rejected`, whose cases
  include the measured inversion (a ploidy-2 set holding one entry) and a haploid keyed
  as a diploid.

### M2 — the two wrong measurements (Applied; upstream deferred)

"29% covered by runs" corrected to 26%, with the realised 0.2629 stated beside the
fitted 0.2634 and the research note's §3.4 cited, so the two numbers can be checked
against each other. `MIN_SITES_TO_FIT`'s precision clause replaced with what the
measurement actually says and what it does not: six million observations is two million
sites at three reads, 200 times this floor, and **what a fit at 10,000 sites is worth
was not measured**.

**`spec/parameter_prepass_generic.md:785` and `arch/parameter_prepass_generic.md:1128`
carry the same 29% slip and are not corrected here** — see §5.

### M3 — the report's own error (Applied)

The stated reason was that a `>=` would let a haploid answer. It would not: `1 >= 2` is
false. The doc comment on the test now says what is true — `>= 2` lets a *tetraploid*
answer and the tetraploid test catches it; `<= 2` lets a haploid answer and only this
test does. The implementation report is corrected, and the commit message of `b9ef1e8`
carries the wrong claim and cannot be corrected without rewriting history.

### M4, M5, Mi7, Mi8 — the untested surfaces (Applied)

Four tests, each verified to fail against the mutation that motivated it:

| mutation | test that fails |
|---|---|
| `Display for Ploidy` body → `write!(f, "MUTANT")` | `ploidy_displays_as_the_bare_copy_number`, plus both message tests |
| `RunsModelStarts::default()`'s two fields swapped | `the_default_starts_are_the_values_the_design_specifies_on_each_axis` |
| length check dropped from `try_new` | `a_set_with_the_wrong_number_of_entries_for_its_ploidy_is_rejected` |

`RunsModelFit::best_start()` added as the single reader of the best-first ordering. The
floor assertions now compare against `MIN_SITES_TO_FIT.to_string()` rather than a
literal, so a floor raised tenfold no longer leaves them green.

### Mi1–Mi6, Mi9, Mi10, N1, N2 — documentation and shape (Applied)

`separations`' doc now leads with the inversion — **a smaller number means the two
states are guessed further apart** — and says that getting it backwards is how a start
set ends up spanning nothing. The `Default` impl is documented and both fields state
their default values. Each floor gained its citation and its soft/fixed label;
`MAX_COUPLED_FIT_ITERATIONS` now justifies itself by what the loop needs (one iteration
on a single-library sample; convergence measured in all 25 worlds) rather than by how
many samples iterate. `MIN_WINDOWS_TO_FIT_INBREEDING` moved to `runs.rs`, which is where
the architecture declares it. The error variants carry `floor`, which deletes the
parent→child import and makes the message right when the STR path raises one.
`ScanResult::log_likelihood` became a `LogProb`. `homozygous_non_reference_rate`'s doc
now says the pair-belonging quantity is heterozygosity *plus* this rate.
`inside_het_floor`'s doc says "floor" names its role rather than its arithmetic.
`by_alt_copies` widened to `SmallVec<[GenotypeFrequency; 5]>`, matching its own doc's
ploidy-4 case. `#[non_exhaustive]`'s justification reworded.

### Deferred

- **M6 — splitting `ParameterEstimationError::Domain`.** Real, and the shape is
  inherited verbatim from `arch` §5.4. Fixing it amends the architecture. The variant's
  doc now states the gap explicitly rather than leaving it implied.
- **M2 upstream — the 29% slip in the spec and the architecture.** Same reason.
- **N3, N4** — `implied_f`'s name and the ten remaining raw `f64`s in `runs.rs`, both of
  which follow the architecture's own declarations.

## 5. Deferred findings to carry forward
M6, the upstream 29% slip, N3 and N4 — all recorded in `PROJECT_STATUS.md`.

## 6. Disputed findings to return to reviewer
None. Two reviewers proposed incompatible fixes for M1; the report records which was
taken and why, which is a resolution rather than a dispute.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code. `SampleRates` gained a validating
constructor, but it is called once per ploidy per sample, not per site.

## 10. Commands run
`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --lib ng::parameter_estimation`, `cargo test --lib ng::`,
`cargo test --all-targets --all-features`, `cargo doc --no-deps --lib`, plus the three
mutations tabled above — all via `./scripts/dev.sh`.

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, no output
- `cargo test --lib ng::parameter_estimation` → 0, 40 passed
- `cargo test --all-targets --all-features` → 101, 2,941 / 1 (pre-existing) / 5
- `cargo doc --no-deps --lib` → 101, 12 unresolved links, all pre-existing
- each of the three mutations → the intended test fails, and only it

## 12. Notes

- **The commit's one recorded assumption was half right, and the wrong half was the one
  that mattered.** Deferring the *sum* check was defensible — nothing reads the sum
  until Milestone E. Deferring the *length* check was not, because an accessor reads by
  dosage today. Two reviewers found it independently, one by probe output and one by
  pricing the fix.
- **Three of the six Majors are wrong numbers rather than wrong code**, and two of those
  are in prose this commit wrote. On a milestone whose deliverable is declarations, that
  is where the defects live.
- **One Major is in the review's own subject's report**, and the correction ("`>=` would
  let a haploid answer") is arithmetic anyone could have checked. It survived writing,
  a commit message, and a chat summary before an agent ran the mutation.
