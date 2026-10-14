# Fix Application Report: str_bundles_reads_psp_2026-10-11.md

**Date:** 2026-10-11
**Source review:** `doc/devel/reports/reviews/str_bundles_reads_psp_2026-10-11.md`
**Source state reviewed against:** `str-bundles` at 3db4e2fc plus its uncommitted diff (step 3)
**Execution mode:** non-interactive (delegated). Two fixes were decided by the owner before the
review (M1, M2); no other question was put to the owner.
**Overall status:** Completed. All seven Majors are applied; the Minors not applied are listed with
their reasons.

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 7
- Minors: 11
- Nits: 1 group

### Outcome totals
- Applied: 9 (M1, M2, M3, M4, M5, M6, M7, Mi1, Mi2)
- Applied with adaptation: 3 (Mi3, Mi4, Mi5)
- Already fixed: 0
- Deferred: 5 (Mi6, Mi8, Mi9, Mi10, Mi11)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 1 (Mi7, absorbed by M1's edit of the same lines)
- Awaiting user answer: 0 (open questions 1 and 2 of the review are left to the owner)

### The two decided fixes, measured

**M1: runs of every motif.** HG002's 2,410 bundles at the calling floors, from
`examples/ng_str_bundle_runs.rs`:

| runs in the bundle | 0 | 1 | 2 | 3 | 4 | 5 | runs in all |
|---|---|---|---|---|---|---|---|
| members' motifs only (before) | 141 | 922 | 1,210 | 124 | 12 | 1 | 3,767 |
| every motif of period 1–6 (after) | 130 | 869 | 1,234 | 155 | 20 | 2 | 3,892 |

The 125 more runs are 9 of period 1, 97 of period 2, and 4, 9, 5 and 1 of periods 3 to 6. The
spec's example, chr1:8,007,443, now comes out `CT:0:40, AT:40:16, T:56:8`, which is
`(CT)20 (AT)8 T8`. Before, it was `CT:0:40, T:55:9`. chr1:28,487,365 is `A:0:9` in both.

**M2: the anchor test.** Per 1,000 crossing reads on the stuttering benchmark (PCR, 30×), measured
by `tmp/step3/spelling.py` as the implementer measured it. The columns are:

- **true:** spelled as one of the sample's two true alleles;
- **stutter shape:** one whole-unit indel away from a true allele;
- **same length:** the same length as a true allele, with another base;
- **part of hull:** spelled over part of the hull only;
- **other length:** any other length.

| | crossing reads | true | stutter shape | same length | part of hull | other length |
|---|---|---|---|---|---|---|
| bundles, 1 sample, before | 10,859 | 921 | 23 | 45 | **7.4** | 3.2 |
| bundles, 1 sample, after | 9,910 | 939 | 16 | 44 | **0.2** | 0.7 |
| segments, 1 sample, after (unchanged) | 59,653 | 932 | 28 | 39 | 0.0 | 0.8 |
| bundles, 50 samples, before | 542,911 | 906 | 38 | 45 | **7.7** | 3.6 |
| bundles, 50 samples, after | 497,333 | 922 | 32 | 44 | **0.4** | 1.7 |
| segments, 50 samples, after (unchanged) | 2,970,633 | 902 | 59 | 38 | 0.2 | 1.6 |

At one sample, 949 reads in 10,859 that were complete observations are now dropped as ending
inside the hull. Reconstructed from the rates (rounded): about 106 were part-of-hull or
other-length spellings, 94 stutter-shaped, 52 the right length with another base, and 693 spelled
as a true allele. The rule demotes a read that anchored on fewer than five flank bases whatever
its spelling. A bundle cell now
has 20 crossing reads and 23 that yield no observation, where it had 22 and 21.

### Validation summary
- `cargo fmt --check` → exit 0 (after `cargo fmt` on the new code).
- `cargo clippy --release --all-targets` → exit 0. No warning in a touched file except the
  pre-existing `type_complexity` at `src/locus_generation/ssr.rs:2448` (it was line 2422 before
  the new tests). One new warning, `manual_is_multiple_of` in `reference_runs.rs`, was fixed.
  `-D warnings` still fails on that line and on `src/run/edge_deletion_gathering.rs:201`, both
  pre-existing.
- `cargo test --release --no-fail-fast --lib -- locus_generation region_typing psp:: run:: calling:: window_coverage parameter_estimation cli:: repeat_catalog alignment::`
  → 4,119 passed, 2 failed, 11 ignored. The 2 failures are the known debug-only `should_panic`
  tests (`alignment::left_align_{repeated,structured}::tests::an_offset_past_the_reference_panics`);
  the third known one is outside this filter. The four `cli::cross_platform_digests` tests pass,
  so the digests did not move.
- `cargo test --release --test ng_repeat_catalog_differential` → 13 passed.
- `cargo doc`, `cargo audit` → not run (no dependency changed).
- Performance check → not applicable: no bench in `benches/` reaches the bundle generator. The
  anchor rule is off for every existing user of the flat-gap aligner.

**End to end on the stuttering benchmark (one sample, PCR, 30×, the fixed build):**

- every non-bundle psp record is byte-identical with `tmp/bins/pvc_main`'s. The md5 of the
  4,120,085 dump rows, `9c81a7f68526eddb9a282f5ee3990968`, is main's.
- the VCF is identical to main's in every column but INFO. INFO differs only in `PARALOG_LR`,
  now in 66 of 2,396 records (67 before the fixes), and the cut moved from 31.75 to 31.85 as
  before.
- the generated psps (`tmp/review_2026-10-11_str-bundles-reads-psp/bench/*/psps`, 731 MB at 50
  samples) were deleted afterwards. The binary is kept as `tmp/bins/pvc_review_fixed`.

### Unresolved high-priority findings
- None of Major severity. Open question 1 of the review (refuse 1.x psps up front?) and question 2
  (should reads ending inside a bundle count toward window coverage? Mi10) are the owner's.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | Runs looked for only in the members' motifs | Apply (decided) | Applied | Owner, before review | `src/region_typing/reference_runs.rs`, `src/locus_generation/ssr_bundle.rs`, `examples/ng_str_bundle_runs.rs`, `str_bundle_loci.md` | Pass | No |
| M2 | Major | Flat-gap delimiter anchors from shape | Apply (decided) | Applied | Owner, before review | `src/alignment/ssr_best_path_flat_gap.rs`, `src/alignment/ssr_unit_robust.rs`, `src/locus_generation/ssr_bundle.rs`, `str_bundle_loci.md` | Pass | No |
| M3 | Major | psp reader accepts runs past the hull | Apply | Applied | No | `src/psp/record.rs` | Pass | No |
| M4 | Major | 2.0 compatibility story stated wrongly | Apply (docs) | Applied | No | `src/psp/header.rs`, `psp_file_format.md`, `str_bundle_loci.md` | Pass | Open question 1 |
| M5 | Major | `--help` says bundles are "not built yet" | Apply | Applied | No | `src/cli/command_line.rs`, `src/cli/call_from_alignments.rs`, `src/run/callers.rs` | Pass | No |
| M6 | Major | Cut-back arithmetic untested | Apply | Applied | No | `src/region_typing/reference_runs.rs` | Pass, 2 mutants killed | No |
| M7 | Major | Generator accounting test vacuous | Apply | Applied | No | `src/locus_generation/ssr_bundle.rs` | Pass | No |
| Mi1 | Minor | `BODY_FIELDS` doc describes the draft layout | Apply | Applied | No | `src/psp/record.rs` | Pass | No |
| Mi2 | Minor | Comments say bundle records are never emitted | Apply | Applied | No | `src/window_coverage/depth.rs`, `src/locus_generation/mod.rs` | Pass | No |
| Mi3 | Minor | Refusal test accepts `Truncated` everywhere | Apply | Applied with adaptation | No | `src/psp/record.rs` | Pass | No |
| Mi4 | Minor | No payload-on-tract refusal test | Apply | Applied with adaptation | No | `src/psp/record.rs` | Pass | No |
| Mi5 | Minor | Right-edge escape untested | Apply | Applied with adaptation | No | `src/locus_generation/ssr.rs` | Pass | No |
| Mi6 | Minor | Flank-only test misses the ordering | Apply | Deferred | No | None | N/A | Yes |
| Mi7 | Minor | `unknown_edges` misnamed | Apply | Superseded by M1 | No | `src/locus_generation/ssr_bundle.rs` | Pass | No |
| Mi8 | Minor | inspect-psp `COLUMNS_VERSION` | Ask | Deferred | No | None | N/A | Yes |
| Mi9 | Minor | `UnsupportedVersion` reads "too new" | Defer | Deferred | No | None | N/A | Yes |
| Mi10 | Minor | Bundle hull reports about half its read depth to coverage | Ask | Deferred | No | None | N/A | Open question 2 |
| Mi11 | Minor | Naming and structure points | Defer | Deferred | No | None | N/A | Yes |
| Nits | Nit | Tie-break, imports, tuple, `by_ref` | — | Partly applied | No | `reference_runs.rs` | Pass | No |

## 3. Questions asked and answers

None asked in this run. The two decided fixes (M1, M2) came with the task. The review's open
questions 1 and 2 are recorded there for the owner.

## 4. Per-finding log

### M1 — Runs looked for only in the members' motifs
- **Severity:** Major
- **Initial decision:** Apply (decided by the owner before the review)
- **Final status:** Applied
- **Reasoning:** the spec's own example lost a 16-base `AT` run to fixed sequence.
- **Implementation summary:** `find_reference_runs(hull, min_copies)` no longer takes the
  members. For every period from 1 to 6, `periodic_stretches_of(hull, period)` yields the maximal
  stretches whose unit passes `is_a_motif_of_its_own_period`. That test requires `A`, `C`, `G` or
  `T` only, so no run of `N` is found, and it rejects a unit that repeats a shorter one, so
  `ATATAT` is found at period 2 only. The longest-first choice, the copy floor and the
  give-up-shared-bases rule are unchanged. The tie-break on the motif's bytes is dropped, because
  it could not decide anything (the review's Nit). `is_rotation_of` is gone: one scan per period
  finds each motif once, in the phase its first base gives it. The bundle generator still takes
  the members for the delimiter's motif. Spec §4, §5.1, §6.1, §10 step 3 and §11 carry
  "Changed in review" notes.
- **Verification performed:** the spec's two §5.1 examples and the reviewer's hand-traced
  fixtures. The HG002 re-run is shown above.
- **Files changed:** `src/region_typing/reference_runs.rs`, `src/locus_generation/ssr_bundle.rs`,
  `examples/ng_str_bundle_runs.rs`, `doc/devel/ng/spec/str_bundle_loci.md`
- **Tests added or modified:**
  - `the_ct_at_t_bundle_comes_out_as_three_runs_…`: now with no `AT` member, which replaces
    `a_motif_no_member_carries_is_not_looked_for`;
  - new: `a_stretch_is_counted_under_its_smallest_period`, `no_run_is_found_in_n_bases`,
    `a_hull_too_short_to_repeat_has_no_run`, `a_unit_is_a_motif_only_at_its_own_period_and_of_acgt`;
  - removed: `rotations_among_the_members_are_one_motif`,
    `rotations_are_recognised_and_other_units_are_not`;
  - the rest now call the finder without members.
- **Validation:** `cargo test --release --lib -- reference_runs` → 17 passed.
- **Residual risk:** the scan can now find a run in a bundle's gap sequence that no member names.
  That is the purpose of the fix, and the copy floor bounds it.

### M2 — The flat-gap delimiter calls a flank anchored from the alignment's shape
- **Severity:** Major
- **Initial decision:** Apply (decided)
- **Final status:** Applied
- **Implementation summary:** `SsrFlatGapAligner` gains `anchor: Option<AnchorRule>`.
  `SsrFlatGapAligner::new` sets `None`, so the parity oracle, the bake-offs and the examples run
  exactly as before. `with_anchor_rule` sets `Some`. With a rule, the traceback counts each
  flank's support and its agreeing bases in the 5-base window beside the junction. It then
  applies the unit-robust `resolve_anchors`, which only demotes and never removes the last anchor.
  A demoted side falls back to the read's own end. `Aligned`, `FlankEvidence` and
  `resolve_anchors` in `ssr_unit_robust.rs` became `pub(super)` so that they are reused rather than
  copied. The bundle generator builds `with_anchor_rule(…, AnchorRule::default())`. Segments use
  the unit-robust delimiter and are untouched.
- **Where the failure came from.** My first fixtures were clean reads, low-quality tails and four
  real hulls, and none made the positional test fabricate a complete span. The benchmark's own
  reads showed the shape: at chr1:144,818 (hull `A10 CTACTC A8`) a read carrying `A11` ends inside
  the hull. The path takes ten `A`s, deletes the rest of the hull, puts the eleventh `A` against
  the right flank's first base and runs off the frame's end. The positional test calls that flank
  anchored on one mismatched base. The unit test uses this read.
- **Review suggestion used verbatim?:** No (the place was chosen as decided: opt-in in the
  aligner).
- **Tests added:** `the_anchor_rule_turns_a_fabricated_complete_into_a_lower_bound`. It asserts
  that the positional test reaches the defect (both flanks anchored, spelled `A10`) and that the
  rule gives a lower bound running to the read's end.
  `the_anchor_rule_leaves_reads_that_hold_their_flanks_unchanged` covers six reads: clean,
  longer, shorter, interrupted, and run off either end.
- **Validation:** the flat-gap module's tests pass. The benchmark table is above.
- **Residual risk:** the rule also demotes reads that were spelled right: at one sample, about
  693 of the 949 newly dropped reads. They anchored on fewer than five flank bases, and they are
  now counted under `ended_inside_a_bundle`. In all, 949 of 10,859 crossing reads at a bundle, about 1 in 11, are now dropped,
  which matters most at low depth (spec §8).

### M3 — The psp reader accepts bundle runs past the hull
- **Severity:** Major
- **Final status:** Applied
- **Implementation summary:** `read_bundle_payload(payload, bytes_in, hull_len)` refuses a run
  whose end passes `hull_len`, the record's `reference_bases.len()`, naming
  `str-bundle-payload`: "a run at {offset} of {length} bases ends past the {hull_len}-base hull".
  Three fixtures had runs past their hulls and now fit them: the every-kind round trip, the
  byte-pinned bundle record and the many-records oracle.
- **Tests added:** `decode_record_body_refuses_a_bundle_run_past_the_hull` (a 7-base hull with a
  run at 5 of 4 bases), and a "past the hull's end" case in the payload refusal table.
- **Validation:** `psp::record` tests pass.

### M4 — The 2.0 compatibility story stated wrongly
- **Severity:** Major
- **Final status:** Applied (documentation; behaviour unchanged)
- **Implementation summary:** three places now say the same thing: `FORMAT_VERSION`'s and
  `OLDEST_READABLE_MAJOR`'s docs, `psp_file_format.md`, and `str_bundle_loci.md` §6.2. They give
  the real reason for 2.0: a 1.x reader would skip the trailing payload and read every bundle with
  no runs, and it would ignore the bundle cap. They also say that a 1.x header passes the version
  check but that every psp written by main must be written again, because it lacks
  `max-str-bundle-len-bp`.
- **Follow-up:** open question 1. Refusing 1.x outright would be a behaviour change and is left to
  the owner.

### M5 — `--help` and docs still say bundles are "not built yet"
- **Final status:** Applied. The `call-from-alignments` help, its module doc and the
  `bundles_set_aside` doc now say that an STR bundle is measured but not yet called, and is
  counted as set aside.

### M6 — Cut-back arithmetic untested
- **Final status:** Applied
- **Tests added:** `a_stretch_cut_at_its_end_goes_back_to_whole_copies` (`T (CAG)8 G30 T` gives
  `(CAG)7` and 21 bases) and `a_stretch_cut_only_at_its_start_keeps_its_part_copy`
  (`T G30 (CAG)8 CA T` gives 26 bases). The old test is renamed
  `a_stretch_gives_up_the_bases_it_shares_with_a_longer_run` for what its fixture shows.
- **Mutation check** (`tmp/review_2026-10-11_str-bundles-reads-psp/mut.sh`):
  - "never cut back" fails 1 test (`mut1.log`);
  - "always cut back" fails 4 tests (`mut2.log`).
  - Before this change both mutants survived. The source was restored and compared with `cmp`.

### M7 — The generator-level accounting test was vacuous
- **Final status:** Applied
- **Implementation summary:** the all-`A` test is replaced by
  `a_bundle_drops_and_counts_the_reads_ending_inside_and_keeps_the_crossing_ones`. Its `chr1` is
  built with a real bundle, `(CA)6 TT (GAT)6` at 41..=72, between flanks that continue neither
  run. Five reads copied from it:
  - two cross the hull;
  - one ends inside the hull;
  - one starts inside the hull;
  - one lies in the left flank alone.
- **What it asserts:**
  - `ended_inside_a_bundle == 2`, `outside_tract == 1`, `observations_partial == 0`;
  - one complete observation of the whole hull, with `num_obs == 2`;
  - the runs `(CA)6` at 0 and `TGA` × 19 bases at 13.
- **Validation:** passes.

### Mi1, Mi2 — Stale docs
- **Final status:** Applied.
  - `BODY_FIELDS`' doc now describes the 24 fields and the two trailing byte strings.
  - `depth.rs` says a bundle record covers its hull, and that its depth counts crossing reads only.
  - "Only one slot is filled today" is gone from `locus_generation/mod.rs`.

### Mi3 — Refusal test accepted `Truncated` for every case
- **Final status:** Applied with adaptation. Each case asserts `Malformed` on
  `str-bundle-payload` with its own reason substring. No deliberate truncation case was added:
  with a run's floor of 4 bytes, a short count is refused by `read_count`, so a truncation case
  would only re-test that.

### Mi4 — Payload on a tract record
- **Final status:** Applied with adaptation. The refusal test now runs over `Generic` and `Ssr`
  records, and its tail check reads the last two bytes, because a tract's body does not end in
  three zeros.

### Mi5 — Right-edge escape
- **Final status:** Applied with adaptation:
  `the_edge_escape_check_asks_the_right_edge_on_its_own` asserts `(None, Some(3))` → escapes and
  `(Some(2), None)` → does not.

### Mi6 — Flank-only test misses the ordering in `partial`
- **Final status:** Deferred
- **Reasoning:** the suggested read (`8M6S`, clip `GGGGGG`) does not do what the finding needs.
  With the clip facing the hull, the delimiter may well place the six clipped bases inside the hull
  as mismatches, and the read is then rightly "ended inside", not "outside". No fixture was found
  that reaches the delimiter and covers no hull position. The behaviour is unchanged.

### Mi7 — `unknown_edges` misnamed
- **Final status:** Superseded by M1, which rewrote these lines. The variable is now
  `edges_the_fetch_ignores`, with a comment saying the fetch reads only the span.

### Mi8, Mi9, Mi11 — inspect-psp version; `UnsupportedVersion` wording; naming and structure
- **Final status:** Deferred.
  - Mi8 is a format-contract decision for the owner: bump `COLUMNS_VERSION` or add a column.
  - Mi9 changes a public error variant.
  - Mi11 covers renames (`tract_len`, `hulls`), visibility, moving `ReferenceRun` out of
    `region_typing`, a separate `RunError` variant, the `Default` placeholder, the `N` fallback,
    payload error offsets and the `FlankClamped` boundary test. Each is a refactor wider than the
    finding, best done with step 4's changes to the same types.

### Mi10 — A bundle hull reports about half its read depth to coverage
- **Final status:** Deferred (open question 2). The coverage at a hull is correct counting, not
  double counting: a hull is one typed region whose positions no other record reports, and before
  this change they were absent from the windows. It is low by the reads that end inside the hull:
  at one sample, 20 crossing reads per bundle cell against 23 dropped.

## 5. Deferred findings to carry forward
- Mi6 — no fixture found that reaches the ordering; the behaviour is unchanged.
- Mi8 — owner's call on the inspect-psp contract.
- Mi9 — a public error variant change.
- Mi10 — owner's call (open question 2).
- Mi11 — refactors to do with step 4.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no. No bench in `benches/` reaches the bundle generator. The flat-gap aligner's
  only new work is behind `anchor.is_some()`, which is false for every existing caller. For a
  bundle locus it adds per read one comparison and two counters per matched flank column of the
  traceback, a walk the aligner already does.

## 10. Commands run
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-reads-psp/{baseline,v1,v3,v4,v5,v6}.sh`
  (fmt, clippy, tests, differential)
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-reads-psp/build.sh`
- `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-reads-psp/hg002_runs.sh`
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-reads-psp/{bench_1s,bench_50s}.sh`
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-reads-psp/mut.sh`

## 11. Command results
- v4 lib tests → exit 101: 4,119 passed, and 2 failed (the known debug-only `should_panic` tests).
- differential → exit 0, 13 passed.
- clippy → exit 0; no new warning in touched files.
- fmt → exit 0.
- the HG002 runs and the spelling tables are in section 1.
- the mutants: 1 and 4 tests failed, so both mutants were killed.

## 12. Notes
- `PROJECT_STATUS.md` was not updated. It has no block for the STR bundle work, as the two earlier
  bundle reviews' reports record.
- Nothing was committed.
