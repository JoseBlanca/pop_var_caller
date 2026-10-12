# Fix Application Report: str_bundles_partial_reads_2026-10-12.md

**Date:** 2026-10-12
**Source review:** `doc/devel/reports/reviews/str_bundles_partial_reads_2026-10-12.md`
**Source state reviewed against:** branch `str-bundles`, HEAD 1ef679aa plus the uncommitted step 4b diff (frozen as `tmp/review_2026-10-12_str-bundles-partials/step4b_under_review.patch`)
**Execution mode:** non-interactive (a subagent of the session lead; owner decisions recorded as Deferred)
**Overall status:** Partial — every Apply finding applied; the two modelling questions and five Minors deferred

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 8 (H1, H2, M2–M7)
- Performance "Likely": 2 (H3, L1)
- Minors: 7 (Mi1–Mi7)
- Nits: 4

### Outcome totals
- Applied: 8 (H1, H2, H3, L1, M5, M6, Mi1, Mi6), and two nits (0.0024; the test wording)
- Applied with adaptation: 1 (Mi7, `lay_out`'s buffers only)
- Already fixed: 0
- Deferred: 8 (M2, M3, M4, M7, Mi2, Mi3, Mi4, Mi5), and two nits (the table header; the file size)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 2 (M2 and M3, recorded as Deferred with the question)

### Validation summary
- `rustfmt --edition 2024 --check` on the eight touched source files → exit 0.
- `cargo clippy --release --all-targets` → exit 0; on the touched files only the pre-existing
  `locus_generation/ssr.rs:2701` `type_complexity` warning remains.
- `cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run::callers window_coverage:: cli::cross_platform_digests`
  → exit 0, `1774 passed; 0 failed; 4 ignored` (includes the four cross-platform digest tests and the
  one-run regression `a_one_run_bundle_scores_as_the_segment_it_is`).
- `cargo test --all-targets`, `cargo doc`, `cargo audit` → not run (the brief asked for touched tests).
- Performance check → not criterion; the cost benchmark of §9 below, with byte-identical calls.

### Unresolved high-priority findings
- M2 — owner decision: should a partial read be neutral about runs it never reached?
- M3 — owner decision: carry the anchored flank so whole-hull reads are not scored at half weight.
- M4 — segments compare right-flank partial reads out of phase (27 in 100 of segment partials); a separate change.
- M7 — no driver-level test that a bundle's partial reads reach calling.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| H1 | Major | Layout rebuilt per read per candidate | Apply | Applied | No | `likelihood/ssr_emission.rs`, `likelihood/ssr.rs` | Pass; VCFs byte-identical | No |
| H2 | Major | Two logarithms per comparison | Apply | Applied | No | `likelihood/ssr_emission.rs` | Pass; VCFs byte-identical | No |
| H3 | Likely | Fixed pieces not memoised | Apply | Applied | No | `likelihood/ssr_emission.rs` | Pass; VCFs byte-identical | No |
| L1 | Likely | Agreeing bases re-added | Apply | Applied | No | `likelihood/ssr_emission.rs` | Pass; bitwise test | No |
| M2 | Major | Unreached runs not neutral | Ask | Deferred | Owner | spec §7.6 note | N/A | Yes |
| M3 | Major | Whole-hull reads at half weight | Ask | Deferred | Owner | spec §7.6 note | N/A | Yes |
| M4 | Major (pre-existing) | Segment right-flank partials out of phase | Defer | Deferred | No | spec §7.6 note | Measured | Yes |
| M5 | Major | Two indel-term tests could not fail | Apply | Applied | No | `likelihood/ssr_emission.rs` (tests) | Pass | No |
| M6 | Major | Partial-read test gaps | Apply | Applied | No | `likelihood/ssr_emission.rs` (tests) | Pass; 2 mutations caught | No |
| M7 | Major | No test that partials reach calling | Defer | Deferred | No | None | N/A | Yes |
| Mi1 | Minor | Debug-only precondition | Apply | Applied | No | `likelihood/ssr_emission.rs` | Pass | No |
| Mi2 | Minor | `HeldFlank::Either` reaches Left/Right helpers | Defer | Deferred | No | None | N/A | With M3 |
| Mi3 | Minor | `reference_hull_len` a bare `u32` | Defer | Deferred | No | None | N/A | Yes |
| Mi4 | Minor | Empty flank on a long deletion | Defer | Deferred | No | None | N/A | Yes |
| Mi5 | Minor | `unreachable!` in `fill_emissions` | Defer | Deferred | No | None | N/A | Yes |
| Mi6 | Minor | Stale documents | Apply | Applied | No | `locus_generation/ssr.rs`, `allele_candidates/ssr.rs`, `str_bundle_loci.md`, `locus_generation_ssr.md` | fmt/clippy pass | No |
| Mi7 | Minor | `lay_out` long; buffer tuples; bound by hand | Apply (part) | Applied with adaptation | No | `likelihood/ssr_emission.rs` | Pass | Rest deferred |
| Nits | Nit | 0.0023; "two bases"; table header; file size | Apply (two) | Applied (two), Deferred (two) | No | spec, test doc | N/A | No |

## 3. Questions asked and answers

1. **M2** — Should a read that ended inside the hull contribute a factor of one for every run it
   never reached (so it is neutral about their counts), or keep each run's total slip probability
   and record the measured 6-in-100 spread in §7.6?
   - **Answer:** awaiting the owner. §7.6 now records the measurement.
2. **M3** — Carry the anchored flank from the read classifier to the emission, so a read covering
   the whole hull is scored from the side it held rather than averaged at half weight?
   - **Answer:** awaiting the owner.

## 4. Per-finding log

### H1 — Layout rebuilt per read per candidate
- **Severity:** Major. **Initial decision:** Apply. **Final status:** Applied.
- **Reasoning:** the profile named the layout and its allocations; none of it depends on the read.
- **Implementation summary:** `CensoredScratch` keeps up to 16 `CensoredLayout`s, each with a
  `CensoredLayoutInputs` compared in full (allele bases, flank beyond the hull, runs, runs' places and
  counts, run models, held flank) before reuse; the layout now also holds the part-repeat renderings
  with their flank (`PartRepeatChange`, `PartRepeatPlace`), built through `part_repeat_places`, the
  place walker `part_repeat_emission` now also uses. The DP is two free functions over the kept
  layout. `fill_emissions` scores partial reads candidate by candidate.
- **Review suggestion used verbatim?:** N/A (the reviewer was the lead).
- **Verification performed:** the touched tests; VCF bodies byte-identical on four runs (§9).
- **Files changed:** `src/calling/likelihood/ssr_emission.rs`, `src/calling/likelihood/ssr.rs`.
- **Residual risk:** a row with more than 16 (candidate, flank, library) combinations lays some out
  again — slower, never different.

### H2 — Two logarithms per comparison
- **Final status:** Applied. `substitution_probability` is now `letters_under(substitution_scores(…))`;
  the censored path calls `substitution_scores` once per read. The complete path is unchanged.

### H3 — Fixed pieces not memoised
- **Final status:** Applied. The memo holds fixed-piece states after the run states
  (`(runs + piece) × read length + offset`).

### L1 — Agreeing bases re-added
- **Final status:** Applied. `MatchedPrefixLetters` with `agreeing_prefix_len`; pinned bit for bit by
  `a_matched_prefix_gives_letters_under_bit_for_bit` (3,444 comparisons under four score sets).

### M2, M3 — what a partial read says
- **Final status:** Deferred (owner). Measurements recorded in `str_bundle_loci.md` §7.6. No code change:
  either fix changes calls.

### M4 — segment right-flank partials
- **Final status:** Deferred to a separate change, as the brief required. Measured with
  `ng_ssr_loci_dump` on one sample (PCR, 30×): 14,175 of 22,661 right-flank partials out of phase.
  Recorded in §7.6.

### M5 — indel-term tests
- **Final status:** Applied. `a_change_one_run_makes_alone_takes_no_indel_term` now reads the
  rendering a wrong indel term would explain and asserts that term would be over a million times the
  expected value; `a_change_no_single_run_makes_is_also_priced_as_a_sequencing_indel` scores two reads,
  each explained by one half of the sum alone, and asserts the explaining half is over a million times
  the other.

### M6 — partial-read test gaps
- **Final status:** Applied. Three tests and a general by-hand sum (`partial_by_hand`, any number of
  pure runs). Mutations run by the lead, one build: the run memo keyed by offset alone → only
  `a_partial_read_over_three_runs_sums_every_slip_set` fails (the old two-run test passes it); the
  whole-hull average without its right reading → only
  `an_either_read_averages_both_readings_where_both_explain_it` fails. 2 run, 2 caught, 0 changed no
  behaviour.

### M7 — partials reach calling
- **Final status:** Deferred: needs a bundle fixture through `call_one_ssr_bundle_locus`; the §10 table
  is the present guard.

### Mi1 — debug-only precondition
- **Final status:** Applied (`assert_eq!`; one integer comparison per call).

### Mi2–Mi5
- **Final status:** Deferred. Each is unreachable today or an interface change; Mi2 belongs with M3.

### Mi6 — stale documents
- **Final status:** Applied, in the five places listed in the review.

### Mi7 — `lay_out`'s shape
- **Final status:** Applied with adaptation: its four working buffers are one named `LayoutBuffers`
  (clippy's `too_many_arguments` after H1). The split into three functions, the bound restated by hand
  and `part_repeat_emission`'s tuple are deferred.

## 5. Deferred findings to carry forward
- M2, M3 — owner decisions.
- M4 — segment change, with its own measurement.
- M7 — driver-level test.
- Mi2–Mi5, the rest of Mi7 — small interface changes.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check

- **Triggered:** yes — the cost was the brief's required fix. Not criterion: the implementer's cost
  benchmark (single thread, 50 samples, PCR, 30×, 500 bundles), CPU time, builds run back to back in
  one container session.
- **Baseline saved:** yes — the build under review, `tmp/bins/pvc_step4b_prereview`, before any edit.
- **Result:** calling the 500 bundles, CPU seconds: step 4 2.0 and 1.9 (session A), 3.0 and 2.2 (B); under review 14.5, 13.9 (A), 23.8, 14.3 (B); reviewed 7.0, 5.7 (A), 9.4, 6.6 (B). Partial reads add about 8–9 ms a bundle at 50 samples, from about 24 ms.
- **Byte identity:** VCF bodies identical before and after on 1 sample 5× (1,898 records), 1 sample 30× (2,591), 50 samples 30× (6,005), HG002 30× (20,347) and the cost benchmark. Scores unchanged: bundles 99/29/80, 207/1/0, 10,400/13/22 (1 false), HG002 862/150/119 (20 false); segment cells 0 differ from main on all three runs.
- **Outcome:** pass.

## 10. Commands run
- `scripts/dev.sh bash tmp/step4b_review/check.sh`
- `scripts/dev.sh bash tmp/step4b_review/cost.sh <label> <binary> <runs>`
- `scripts/dev.sh bash tmp/step4b_review/perf.sh <binary>`
- `DEV_EXTRA_MOUNT=…/benchmarks scripts/dev.sh bash tmp/step4b_review/bench.sh`
- `scripts/dev.sh bash tmp/step4b_review/dump.sh`; `uv run --no-project python tmp/step4b_review/segment_right_partials.py …`

## 11. Command results
See §1's validation summary and §9.

## 12. Notes
- The final binary is `tmp/bins/pvc_step4b_reviewed`.
- The psps generated for these measurements were deleted afterwards.
