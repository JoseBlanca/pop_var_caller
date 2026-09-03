# Fix Application Report: ng_psp_mode_a1_2026-09-03.md

**Date:** 2026-09-03
**Source review:** `doc/devel/reports/reviews/ng_psp_mode_a1_2026-09-03.md`
**Source state reviewed against:** branch `ng-psp-mode`, uncommitted diff over 9aa05795
**Execution mode:** non-interactive (plan-driven step A1 loop)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 6 (M1–M6)
- Minors: 13 (Mi1–Mi13)
- Nits: 1 grouped set

### Outcome totals
- Applied: 16 (B1, M1, M2, M3, M4, M5, Mi1, Mi2, Mi3, Mi4, Mi5, Mi6, Mi7, Mi8, Mi10, Mi12, Mi13 — Mi6 and Mi13 bundled into M5's rework, see log)
- Applied with adaptation: 1 (Mi11 — one rename adapted)
- Deferred: 2 (M6, Mi9 — both routed to the Checkpoint A conversation, neither blocks)
- Disputed: 0
- Failed validation: 0
- Nits: `# Errors` heading and `to_span` rename applied; `eprintln` kept deliberately (the
  harness captures it on passing runs, and it carries the measured headroom figure); import
  grouping in bench/examples left as-is (the `no_catalog` change removed most of the inline
  paths anyway).

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (one
  `field_reassign_with_default` hit during the run, fixed before completion)
- `cargo test --all-targets --all-features` → 0, all 16 test binaries pass; lib suite
  **6,050 passed; 0 failed; 14 ignored** (was 6,044 pre-fixes; `ng::psp` 415, `region_typing` 96)
- `cargo doc --no-deps` → not run (no public-API doc change beyond doc comments; carried by
  clippy's doc lints)
- `cargo audit` → not run (no dependency change)
- Performance check → **skipped**: every changed line runs once per file at header
  encode/decode or in test/fixture code; the benches time the per-record walk/write loops,
  which no fix touches. No baseline was saved, consistent with the skip.

### Unresolved high-priority findings
- None. (M6 is deferred by design — an owner decision at Checkpoint A, not an open defect.)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | fixture catalog sub-sections are pure defaults | Apply | Applied | segmentation_section.rs | Pass |
| M1 | Major | newline in catalog strings forges header lines | Apply | Applied | segmentation_section.rs, header.rs | Pass |
| M2 | Major | empty-analysed-regions rule untested | Apply | Applied | header.rs | Pass |
| M3 | Major | u64→u32 bound clamp untested | Apply | Applied | segmentation_section.rs | Pass |
| M4 | Major | whole-contig span + TOML-ceiling untested | Apply | Applied | segmentation_section.rs, header.rs | Pass |
| M5 | Major | diff edits frozen src/regions.rs | Apply | Applied | regions.rs (reverted), region_typing/mod.rs, segmentation_section.rs | Pass |
| M6 | Major | psp↔run dependency direction | Defer | Deferred | None | N/A |
| Mi1 | Minor | criteria bounds unenforced (purity, bundle, period-max) | Apply | Applied | segmentation_section.rs | Pass |
| Mi2 | Minor | wire_float_of comment cites a wrong invariant | Apply | Applied | segmentation_section.rs | Pass |
| Mi3 | Minor | ZeroMin arm untested, is_err-only assertion | Apply | Applied | segmentation_section.rs | Pass |
| Mi4 | Minor | wire keys unpinned against rename | Apply | Applied | header.rs | Pass |
| Mi5 | Minor | "no catalog" sentinel triplicated and unnamed | Apply | Applied | repeat_catalog/mod.rs, bench + 2 examples | Pass |
| Mi6 | Minor | bounds-clamp expression at four sites | Apply | Applied (in-crate sites) | segmentation_section.rs | Pass |
| Mi7 | Minor | for_period(u8::MAX) sentinel probe | Apply | Applied | segment_criteria.rs, segmentation_section.rs | Pass |
| Mi8 | Minor | non-exhaustive encode-side field access | Apply | Applied | header.rs, segmentation_section.rs | Pass |
| Mi9 | Minor | required section in unchanged format 1.0 | Dispute-as-ruled | Deferred (record at checkpoint) | None | N/A |
| Mi10 | Minor | indexed parallel lists in check_catalog | Apply | Applied | segmentation_section.rs | Pass |
| Mi11 | Minor | naming cluster | Apply | Applied with adaptation | region_typing/mod.rs, segmentation_section.rs | Pass |
| Mi12 | Minor | catalog-md5 refusal names no contig | Apply | Applied | segmentation_section.rs | Pass |
| Mi13 | Minor | fixture panics on contigs under 151 bp | Apply | Applied | segmentation_section.rs | Pass |

## 3. Questions asked and answers

None asked mid-run. Two items carry to the Checkpoint A conversation instead of blocking the
step: M6 (where `SegmentationInputs` should live once psp and run are mutually dependent) and
Mi9 (recording the "no version bump — no file predates A1" premise as a deliberate ruling).

## 4. Per-finding log (condensed; one entry per material decision)

### B1 — fixture defaults
Applied as suggested in spirit, adapted in values: `built_under` differs in `min_flank_bp`
(11) and `classification.min_score` (9); `scan` differs in all three fields (3/5/4 against
defaults 2/7/2). A new test `the_fixture_differs_from_every_default_the_operands_have` pins
the claim itself, so the fixture cannot silently regress to defaults again; every round-trip
test inherits the discrimination (this is what moves the lib count 6,044 → 6,050 together
with the other new tests). Test-first: the review's own mutations (M-A/M-B) stand as the
demonstration; the strengthened fixture is the fix.

### M1 — forged header lines
`check_catalog` now runs `check_plain_single_line` over `tool-version` and every catalog
contig name — the same rule, stated for the same reason, as the header's contig-name and
manifest-field-name checks. Two rows added to the two-sided rule table in `header.rs`
(`every_rule_is_refused_by_the_writer_and_by_the_reader_alike`) so writer and reader are both
pinned, plus two direct cases in the section's own rule test.

### M2 / M3 / M4 — the coverage majors
- M2: a rule-table row ("a segmentation recording no analysed ground at all") pins both sides.
- M3: `a_span_on_a_contig_longer_than_u32_max_round_trips` — a contig of `1 << 32`, where the
  truncating cast the mutation introduced gives a zero bound and refuses the writer's own file.
- M4: `a_whole_genome_header_round_trips` (full header path, `header.rs`),
  `a_whole_contig_analysed_span_round_trips` (section level), and
  `a_criteria_value_at_the_toml_ceiling_is_accepted_and_one_more_refused`
  (`min_flank_bp` at exactly `i64::MAX`, then one more).

### M5 — the frozen-file fix
`src/regions.rs` restored byte-identical to `HEAD` (`git restore`). `GenomeRegions` now owns
its span list: `whole_contigs`/`from_bed_path` still parse through `RegionSet`'s public API
and copy the result out; the new `from_normalized_spans` (validation moved from the deleted
`regions.rs` code, semantics unchanged) is wholly ng's. A `spans()` accessor exposes the
`u32` shape the psp encode needs. The three unit tests moved with it and a new property test
(`from_normalized_spans_accepts_every_normalized_set`) pins acceptance against the BED
parser's own normalization as oracle — closing the review's reliability Minor about the
missing property test at the same time. Consumers were verified to use only
`whole_contigs`/`from_bed_path`/`iter`/`len`/`is_empty`, so nothing else changed.

### Mi1 — criteria bounds
`check_criteria` now enforces the classifier's own three release-asserted invariants: purity
in `[0, 1]` (subsumes the NaN arm), `bundle_threshold >= 1`, `periods.max() <= MAX_MOTIF_LEN`.
Each has a field-asserting refusal case in the section's rule test (purity 7.0, NaN, bundle 0,
period-max 9 via the catalog's `built-under`).

### Mi5 — the sentinel named
`RepeatCatalogHeader::no_catalog(tool_version)` with a doc that says what it is (the identity
of no catalog at all, for synthetic files; two such files comparing equal is correct); the
bench and both examples now call it, deleting three 8-line literals.

### Mi11 — naming (adapted)
`from_genomic_order_spans` → `from_normalized_spans` (the reviewer's suggestion);
`contigs: usize` → `contig_count`; `as_span` → `to_span`; the test binding `one_contig` →
`with_a_shortened_second_contig`; the wire doc no longer says "the *walk* asked" (now "this
run asked … the same value `SegmentationInputs::repeat_tract_criteria` holds in memory",
matching `segments.rs`). Adaptation: the reviewer's `previous`/`ordered` nit became
`previous_span`/`in_order` during the M5 move rather than a separate patch.

### Deferred
- **M6** — moving `SegmentationInputs`'s file is an architecture edit (the arch's §4 sketch
  itself puts the type inside the psp header); raised at Checkpoint A with the review's
  proposed `src/ng/segmentation_inputs.rs` + re-export shape as the candidate.
- **Mi9** — no code change is right if the plan's premise stands (no psp outside tests
  predates A1); the premise and its consequence (a pre-A1 *scratch* file refuses with a
  missing-field message, not a version message) are put before the owner at Checkpoint A so
  the decision is recorded rather than accidental.

## 5. Deferred findings to carry forward
- M6 — `SegmentationInputs`'s module home (Checkpoint A).
- Mi9 — record the no-version-bump ruling (Checkpoint A).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — every changed production line is on the once-per-file header
  encode/decode path or in test/fixture/bench-setup code; the bench harnesses time the
  per-record walk and write loops, which no fix touches.
- **Baseline saved:** no (consistent with the skip).

## 10–11. Commands run and results
- `./scripts/dev.sh cargo fmt --check` → 0
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` → 0
- `./scripts/dev.sh cargo test --lib 'ng::psp'` → 0, "415 passed; 0 failed"
- `./scripts/dev.sh cargo test --lib 'region_typing'` → 0, "96 passed; 0 failed"
- `./scripts/dev.sh cargo test --lib 'regions::'` → 0, "61 passed; 0 failed" (production's own
  count — the reverted file's tests are exactly what they were)
- `./scripts/dev.sh cargo test --all-targets --all-features` → 0, 16 binaries, lib
  "6050 passed; 0 failed; 14 ignored"

## 12. Notes
- One mid-run failure worth recording: the new property test's first draft synthesized BED
  spans that could *start* past the shorter contig's end, which the parser rightly refuses
  (`IntervalBeyondContig`); the generator now folds the start into the contig. The failure was
  the generator's, not the code under test.
- The 30,000-scaffold worst-case header re-measures at **10,798,518 bytes** of the
  16,777,187-byte ceiling after the fixture changes (the impl report and PROJECT_STATUS are
  updated to this figure).
