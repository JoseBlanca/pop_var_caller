# Fix Application Report: ng_psp_mode_c1_2026-09-03.md

**Date:** 2026-09-03
**Source review:** `doc/devel/reports/reviews/ng_psp_mode_c1_2026-09-03.md`
**Source state reviewed against:** f00d56e9 + the uncommitted C1 diff, branch `ng-psp-mode`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 4
- Majors: 7
- Minors: 18
- Nits: 2

### Outcome totals
- Applied: 4 Blockers, 6 Majors, 9 Minors, 1 Nit
- Deferred: 1 Major (M7), 9 Minors, 1 Nit — §5
- Applied with adaptation / Already fixed / Disputed / Failed validation / Blocked / Superseded: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib 'pop_var_caller_exp'` → 0, **112 passed** (was 100: 88 direct-mode + 24 here, up from 12)
- **Three of the review's surviving mutations re-run against the fixed tree and each now killed by a named test:** the catalog checked against the digest-free reference view → `a_catalog_built_on_another_reference_of_the_same_shape_is_refused`; a swallowed walk failure → `a_walk_that_stops_names_its_sample_and_leaves_the_earlier_samples_psps_written`; `--min-purity` severed from the routing → `the_min_purity_flag_reaches_the_criteria_the_catalog_is_asked_with`. The file was restored byte-identically after each (verified by `diff`).
- **Re-run for real** after the fixes: the same tomato slice wrote the same `SRS3394712.psp`, **914,715 bytes**, exit 0, with no `.partial` left behind.
- Performance check → skipped: nothing on a `benches/` path changed.

### Unresolved high-priority findings
None. M7 is deferred as a design question for the owner (§5).

## 2. Findings table

| ID | Severity | Title | Final status |
|---|---|---|---|
| B1 | Blocker | psp name built from `@RG SM` unchecked | Applied |
| B2 | Blocker | a stopped re-walk destroys the psp it replaces | Applied |
| B3 | Blocker | catalog-vs-reference digest check untested | Applied |
| B4 | Blocker | mid-cohort failure path untested | Applied |
| M1 | Major | four flags never given a non-default value | Applied |
| M2 | Major | `--min-purity` absent from the cross-command comparison | Applied |
| M3 | Major | timestamp silently falls back to the epoch | Applied |
| M4 | Major | output directory judged last, not first | Applied |
| M5 | Major | `.psp` extension unpinned | Applied |
| M6 | Major | no command-level test walks a sample with reads | Applied |
| M7 | Major | hard-coded filter/generator values invisible at the surface | Deferred |
| Mi1–Mi9 | Minor | see §4 | Applied |
| Mi10–Mi18 | Minor | see §5 | Deferred |

## 3. Questions asked and answers
None — non-interactive. The review's two open questions were resolved as recorded there: a
sample name that is not a file name is **refused at the door**; the hard-coded read-filter and
generator values stay hard-coded (M7, deferred to the owner).

## 4. Per-finding log (applied)

- **B1 — Applied.** `refuse_a_sample_name_that_is_not_a_file_name` requires exactly one normal
  path component and is run for every sample **before the reference is read**; new variant
  `SampleNameNotAFileName`. Pinned by a test over six rejected names (`../elsewhere`, `lane/1`,
  empty, `.`, `..`, `/absolute`) and four accepted.
- **B2 — Applied.** Each sample is walked into `<sample>.psp.partial` and renamed once whole;
  a stopped walk removes the stump and leaves the final path untouched. Pinned by
  `a_stopped_rewalk_does_not_destroy_the_psp_it_was_replacing`, which walks the cohort, breaks
  one sample's input, re-walks it into the same directory and compares the original psp byte
  for byte.
- **B3, B4, M1, M5 — Applied**, using the reliability agent's seven verified test bodies
  essentially as supplied (the catalog test, the stopped walk, `--regions`, `--catalog`,
  `--min-purity`, `--build-index-if-missing`, and the psp names read from the directory rather
  than recomputed).
- **M2 — Applied.** The cross-command test now compares the **whole** `StrRepeatCriteria`
  against what direct mode's own defaults produce, instead of a four-term tuple that could go
  stale as axes are added. That closes the `--min-purity` default hole and the `--min-copies`
  literal duplication together.
- **M3 — Applied.** `provenance()` returns `Result` and refuses an unparseable stamp with a new
  `Timestamp` variant, matching `estimate_contamination.rs` and `cli.rs`.
- **M4 — Applied.** The output directory and the run's timestamp are settled first, then the
  read groups and every sample's name, then the reference, then the ground — and the `# Errors`
  paragraph now names that order.
- **M6 — Applied.** The shared fixture gives zeta three reads and leaves alpha empty, so one
  sample exercises a walk that produces records and the other §12.9's analysed-but-empty case;
  `a_samples_reads_reach_its_psp` asserts both.
- **Mi (applied):** `files_of` → `alignment_files_of`, reading the sample's own read-group id
  list instead of re-scanning the table by string; the subcommand name is one `pub const
  SUBCOMMAND` used by both the header and a test that parses it through clap, so renaming the
  variant breaks the test; `Walk`'s two-state ambiguity **dissolved** by B2's rename scheme
  (the named path never holds a partial file) and its doc now says so; the false present-tense
  claim in `gatherer.rs` corrected to describe what the command actually does; `--regions`'
  forward reference to a command that does not exist removed; `--max-str-len` no longer
  promises a report; `--help` and the module doc now say what a psp *is*, and name the two
  things this stage does not do yet.

## 5. Deferred findings to carry forward

- **M7 — the eleven hard-coded behavioural values** (read filters, five generator knobs) are
  invisible at the command surface. Deferred as a design question: direct mode hard-codes the
  same values, so exposing them is a change to both surfaces and the owner's call. The read
  filters and `max_record_span` are recoverable from the psp header; the other four generator
  knobs are recoverable from nothing, which is the part worth raising.
- **The reference-open block duplicated verbatim with the sibling** (22 lines) — `run_ground`'s
  doc already claims to own "the assembly every mode does before a read is decoded", so this
  belongs there. Deferred to keep C1's diff to C1; it is the natural companion to the C2/C3
  work in the same file.
- **`#[command(flatten)]` over `RepeatRouting`** for the five routing flags — would delete two
  copies of the `#[arg]` blocks, both `ground_request`s and the drift test's reason to exist.
  Deferred at Medium confidence, per the smells agent's own note that no third consumer is
  coming; the live hazard it guards against is now closed by M2's whole-struct comparison plus
  the non-default flag tests.
- **The `a_cohort_on_disk` fixture duplicated** between the two commands' test modules (49
  lines).
- Smaller: the four remaining untested error variants (`Reference`, `ReferenceVerification`,
  `OutputDir`, `Timestamp`); `psp_path_for`'s `pub` with no external caller; the test message
  claiming more than it asserts; the by-construction equality assertion; an empty output
  directory left by a failed run; `GeneratePspsArgs` lacking `Clone`; the four shape-named
  locals inherited from the sibling.

## 6–8. Disputed / Failed validation / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched code reachable from `benches/`.

## 10–11. Commands run and results
In the container via `scripts/dev.sh`: `cargo fmt` / `--check` → 0; `cargo clippy
--all-targets --all-features -- -D warnings` → 0; `cargo test --lib 'pop_var_caller_exp'` → 0
(112 passed, 1 ignored); three mutation re-runs, each failing exactly one named test and the
file restored byte-identically after; and the real run quoted in §1.

## 12. Notes
B2 changes what the command leaves on disk after a failure, which is behaviour C3 will build
on: with the rename scheme in place, C3's overwrite refusal is a check on the final path
before the walk starts, not a rescue after it.
