# Fix Application Report: str_partial_read_phase_2026-10-12.md

**Date:** 2026-10-12
**Source review:** `doc/devel/reports/reviews/str_partial_read_phase_2026-10-12.md`
**Source state reviewed against:** branch `str-partial-phase`, a87ce7b6 plus the uncommitted change
**Execution mode:** non-interactive (a sub-agent; no questions could be asked)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 2
- Minors: 10
- Nits: 4 groups

### Outcome totals
- Applied: 3 (B1, M1, Mi8)
- Applied with adaptation: 2 (M2, Mi7)
- Already fixed: 0
- Deferred: 8 (Mi1–Mi6, Mi9, Mi10)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` (the two touched source files) → 0, clean.
- `cargo clippy --release --all-targets -- -D warnings` → 101, only on two pre-existing lints in
  untouched files (`src/locus_generation/ssr.rs:2117` `type_complexity`,
  `src/run/edge_deletion_gathering.rs:201` `single_match`); nothing in the touched files.
- `cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run::callers
  cli::cross_platform_digests` → 0, `1674 passed; 0 failed; 3 ignored`. The four cross-platform
  digest tests are among them and pass **unchanged**: the fixes move no checksum.
- `cargo doc`, `cargo audit` → not run (no public-doc or dependency change of note).
- Performance check → not applicable (no `benches/` harness reaches the STR emission's right-held
  branch differently; the added work is one rotation check per right-held call).

### Unresolved high-priority findings
- None. M2's fix catches a calling length **longer** than the one the witnesses were laid on, not a
  shorter one (residual risk, below).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | Candidate ending mid-copy mirrored out of phase | Apply | Applied | No | `src/calling/likelihood/ssr_emission.rs`, `doc/devel/ng/spec/read_likelihoods.md` | Pass | No |
| M1 | Major | Symmetric interrupted mirror test cannot fail | Apply | Applied | No | `src/calling/likelihood/ssr_emission.rs` | Pass | No |
| M2 | Major | Neither-border witness silently blended; length unchecked | Apply | Applied with adaptation | No | `src/calling/likelihood/ssr.rs`, `doc/devel/ng/spec/read_likelihoods.md` | Pass | Yes |
| Mi1 | Minor | `of_witness` has no direct test | Defer | Deferred | No | None | N/A | Yes |
| Mi2 | Minor | `HeldFlank` duplicates `AnchoredBorder` | Defer | Deferred | No | None | N/A | Yes |
| Mi3 | Minor | Naming (closure `from`, `witnessed`, `tract_len`, `main` bindings) | Defer | Deferred | No | None | N/A | Yes |
| Mi4 | Minor | Weight one half is an unnamed literal | Defer | Deferred | No | None | N/A | Yes |
| Mi5 | Minor | `fill_emissions` over 120 lines | Defer | Deferred | No | None | N/A | Yes |
| Mi6 | Minor | `TEST_TRACT_LEN` of 1,000 | Defer | Deferred | No | None | N/A | Yes |
| Mi7 | Minor | Digest notes lack the isolating check | Apply | Applied with adaptation | No | `src/cli/cross_platform_digests.rs` | Pass | No |
| Mi8 | Minor | Impl report's homopolymer explanation half wrong | Apply | Applied | No | `doc/devel/reports/implementations/str_partial_read_phase_2026-10-12.md` | N/A (doc) | No |
| Mi9 | Minor | `hg002_two_binaries.sh` hides failures; one build dir | Defer | Deferred | No | None | N/A | Yes |
| Mi10 | Minor | No reused-scratch test across sides | Defer | Deferred | No | None | N/A | Yes |

Nits: not acted on.

## 3. Questions asked and answers

None. The run was non-interactive; nothing needed a policy decision.

## 4. Per-finding log

### B1 — A candidate that ends part-way through a copy is mirrored out of phase
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** reproduced with a failing test before any production edit: at 8 bases the
  right-held read of `CACACAC` scored 3.373e-3 against its mirror image's 3.541e-3.
- **Implementation summary:** new `motif_phased_to_the_end_reversed(candidate, motif, out)`: the
  candidate's last `period` bases, written backwards, when they are the motif in some phase; else
  the motif reversed. The right-held branch of `censored_emission` uses it instead of reversing the
  motif. For a candidate of whole copies — every reference tract, and every pure candidate — the
  result is the motif reversed, as before. Doc comments on `MirroredCensoredRead`,
  `StutterSubstitutionScratch::mirrored` and `censored_from_the_tract_start`, and the spec note in
  `read_likelihoods.md` §5.2, say so.
- **Review suggestion used verbatim?:** Yes, in substance.
- **Adaptation:** None. A rotation by the candidate's length was considered and rejected: it
  splits interrupted tracts whose interruption shifts the phase badly (`CACATCACA` would segment as
  runs of one); reading the phase off the last copy does not.
- **Verification performed:** test fails before, passes after.
- **Files changed:** `src/calling/likelihood/ssr_emission.rs`, `doc/devel/ng/spec/read_likelihoods.md`
- **Tests added or modified:** `a_right_held_read_of_a_candidate_ending_mid_copy_mirrors_the_left_held_score`
- **Validation:**
  - targeted test before the fix → 101, `8 bases: right 0.00337285823656753, mirror from the left 0.003541442065278368`
  - `cargo test … calling::likelihood::ssr_emission` after → 0, `32 passed`
  - the full touched-module run (§1) → 0
- **User input:** None
- **Follow-up:** None
- **Residual risk:** a candidate ending on an interruption falls back to the motif reversed; such a
  candidate has no phase to read, and the part-repeat tiling past it is a convention either way.

### M1 — The interrupted mirror test cannot fail for what it claims
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** folded into B1's test, which adds two interrupted candidates at periods 2 and 3
  that are not palindromes (`CACACAGCACA`, `CAGCAGCATCAGCA`, asserted asymmetric in the test), with
  the mirror image's motif written out by hand (`AC`, `ACG`), and reads that reach past the
  candidate so part-repeat terms are exercised.
- **Implementation summary:** as B1. The old test is kept; its doc comment still overstates it
  (left for the author, Nit-level once the new test exists).
- **Review suggestion used verbatim?:** No — merged into B1's test rather than a separate one.
- **Adaptation:** one test serves both findings.
- **Verification performed:** passes after B1; its first case failed before B1.
- **Files changed:** `src/calling/likelihood/ssr_emission.rs`
- **Tests added or modified:** as B1.
- **Validation:** as B1.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M2 — A neither-border witness is silently scored half from each side
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** the state is unreachable from the generator, so meeting it means the tract length
  and the witness disagree; it must not be scored.
- **Implementation summary:** in `fill_emissions`, an `assert!` (release, not debug-only) that a
  partial whose side is `None` is flush with both borders, with a message naming the cause. The
  spec note says a neither-border witness stops the row.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** the reviewers proposed either a single length source or an assertion; the
  assertion is the narrow fix. A single source is part of the str-bundles merge (one length source
  for segments and bundles), so it is deferred there.
- **Verification performed:** the new `should_panic` test lays a right-held witness on 8 positions
  and scores it against 10; the benchmark runs below ran without a panic.
- **Files changed:** `src/calling/likelihood/ssr.rs`, `doc/devel/ng/spec/read_likelihoods.md`
- **Tests added or modified:** `a_partial_witness_touching_neither_border_is_refused`
- **Validation:** `cargo test … calling::likelihood::ssr` → 0, `58 passed`; full run (§1) → 0.
- **User input:** None
- **Follow-up:** one locus-length source when str-bundles is merged.
- **Residual risk:** a calling length **shorter** than the witness's makes a long left-held read
  look whole-tract (blended) and a long right-held read look left-held. Not caught. Today the two
  lengths are the same quantity (the record region is the segment; pinned for one sample by
  `the_mints_locus_length_is_the_emitted_regions`).

### Mi1 — `of_witness` has no direct test
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** not required by a higher fix; str-bundles replaces the function with
  `of_witness(witness, locus_len, recorded)`, which is where its test belongs.
- **Follow-up:** test the merged function's five answers.

### Mi2 — `HeldFlank` duplicates `AnchoredBorder`
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** touches the generator, which this change does not; do it on the str-bundles
  rebase, where both enums meet.

### Mi3 — Naming
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** renames across the trait and its implementors collide with str-bundles' own
  signature change; settle names once, in the merged signature.

### Mi4 — Weight one half is an unnamed literal
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** the mixture moves into one place on the str-bundles merge; name it there.

### Mi5 — `fill_emissions` over 120 lines
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** str-bundles reorders the same loop; extract after the merge.

### Mi6 — `TEST_TRACT_LEN` of 1,000
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** it hides no existing right-held coverage (checked: no fixture under `src/calling/`
  passing through `SsrSampleEvidence` builds one), and after M2 a future right-flush fixture panics
  rather than blending.

### Mi7 — Digest notes lack the isolating check
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied with adaptation
- **Implementation summary:** one sentence added to each 2026-10-12 note: scoring every partial
  from the left again restores the old checksum, and the review's fixes did not move it. The
  real-cohort measurement the module header asks for is the stuttering benchmark and HG002 of the
  implementation report, cited rather than repeated.
- **Validation:** digest tests pass (§1).

### Mi8 — The homopolymer explanation
- **Severity:** Minor · **Initial decision:** Apply · **Final status:** Applied
- **Implementation summary:** implementation report §4 corrected in place, marked as a review
  correction, with the trace (18 of 20 period-1 verdict changes at 50 samples, 5×, have an
  interrupted true allele).

### Mi9 — `hg002_two_binaries.sh`
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** the script is byte-identical to the copy str-bundles committed in 199541ac; fix it
  once, on main after both land, so the two copies do not diverge.

### Mi10 — No reused-scratch test across sides
- **Severity:** Minor · **Initial decision:** Defer · **Final status:** Deferred
- **Reasoning:** every buffer is cleared before use today; worth a test, not required by a higher
  fix.

## 5. Deferred findings to carry forward
- Mi1, Mi4 — on the str-bundles rebase (one side function, one mixture).
- Mi2, Mi3, Mi5 — on the str-bundles rebase (shared enum, merged signature, reordered loop).
- Mi6 — low value after M2.
- Mi9 — once both branches are on main.
- Mi10 — a test to add.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — the only added work is a rotation check over at most six bases per right-held censored
call; not reached differently by any bench in `benches/`.

## 10. Commands run
- `scripts/dev.sh tmp/review_test_mid_copy.sh ending_mid_copy` (before and after B1)
- `scripts/dev.sh tmp/review_test_mid_copy.sh calling::likelihood::ssr_emission`
- `scripts/dev.sh tmp/review_test_mid_copy.sh calling::likelihood::ssr`
- `scripts/dev.sh tmp/review_verify.sh` (fmt check, the touched-module tests with the digests,
  clippy)
- `scripts/dev.sh tmp/review_clippy.sh` (twice)
- `scripts/dev.sh tmp/review_bench.sh` (the re-measurement below)

## 11. Command results
- see §1 and the per-finding logs.

## 12. Notes

### Re-measurement

B1 changes scores, so two runs of the stuttering benchmark were re-measured: true stutter, the
three binaries calling the same freshly made psps. *main* is a87ce7b6, *implementer* the change as
reviewed, *reviewed* this run's build. Cells as in the implementation report (variant cells right /
wrong / missed, then false calls on reference cells; the three compound hexanucleotide motifs
left out).

| run | main | implementer | reviewed |
|---|---|---|---|
| PCR, 1 sample, 15× | 1,031 / 47 / 47 / 17 | 1,033 / 45 / 47 / 17 | 1,033 / 45 / 47 / 17 |
| PCR, 50 samples, 5× | 40,074 / 8,257 / 8,290 / 1,441 | 40,196 / 8,179 / 8,246 / 1,406 | 40,196 / 8,179 / 8,246 / 1,406 |

The reviewed build's VCF records are identical to the implementer's in both runs (0 differing
lines below the header), and no run panicked on M2's assertion. That is expected: the simulated
reference tracts are trimmed to whole copies, so B1 reaches only read-proposed candidates that end
mid-copy, and on those only the part-repeat terms of right-held reads. The main and implementer
columns reproduce the implementation report's numbers exactly. The copied runs were deleted
afterwards; the implementer's `score_cells.tsv` were put back.
