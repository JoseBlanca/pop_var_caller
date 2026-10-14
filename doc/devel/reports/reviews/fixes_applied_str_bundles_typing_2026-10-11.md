# Fix Application Report: str_bundles_typing_2026-10-11.md

**Date:** 2026-10-11
**Source review:** `doc/devel/reports/reviews/str_bundles_typing_2026-10-11.md`
**Source state reviewed against:** `str-bundles` at 90e6f69e plus the uncommitted step-2 diff
**Execution mode:** non-interactive
**Overall status:** Partial — correctness fixes applied; one owner decision (M2) and refactors deferred

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 6
- Minors: 17
- Nits: grouped

### Outcome totals
- Applied: 10 (M1, M3, M5, M6, Mi3, Mi4, Mi5, Mi6, Mi7, Mi9) and Mi10 with M5
- Applied with adaptation: 1 (M2: documented, behaviour unchanged pending the owner)
- Deferred: 12 (M4, Mi1, Mi2, Mi8, Mi11–Mi17)
- Disputed / Failed / Blocked / Superseded: 0
- Awaiting user answer: M2 and Mi8 (review questions 1 and 2)

### Validation summary (dev container, release)
- `cargo test --release --no-fail-fast --lib` → exit 101: **4,996 passed**, 3 failed — the three
  debug-only `should_panic` tests that fail in release, the same three as before any fix
  (baseline 4,995 passed; one test added to the library here, the others extend existing ones).
- `cargo test --release --test ng_repeat_catalog_differential` → exit 0, **13 passed** (12 before;
  the new one is M1's).
- `cargo test --release --lib -- cli::cross_platform_digests …` → exit 0 after re-recording.
- `cargo fmt --check` on the 14 touched Rust files → clean.
- `cargo clippy --release --all-targets -- -D warnings` → exit 101 on two `type_complexity`
  errors in files this change does not touch (`locus_generation/ssr.rs:2117`,
  `run/edge_deletion_gathering.rs:201`); none in a touched file. The lint stops at the library, so
  the integration-test target was not linted.
- `cargo doc`, `cargo audit` → not run (no public-API docs restructured, no dependency change).
- Performance check → skipped: no bench covers `rows_for`; on real chromosomes the longest tract
  exceeds the cap, so the padded window does not change (HG002's typing is byte-identical).

### Unresolved high-priority findings
- **M2** — overlapping and nested members; owner decision before step 3.
- **M4** — the catalog route's no-clean-cut fallback is still untested.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | windowed read types a chain as a bundle | Apply | Applied | No | `repeat_catalog/mod.rs`, `reader.rs`, differential test | Pass + mutation | No |
| M2 | Major | members overlap and nest | Ask | Applied with adaptation (docs only) | Yes | `region_typing/mod.rs`, `str_bundle_loci.md`, `typed_regions.md`, glossary | N/A | Yes |
| M3 | Major | hull not pinned to detected spans | Apply | Applied | No | `region_typing/mod.rs` (test) | Pass | No |
| M4 | Major | catalog no-cut fallback untested | Apply | Deferred | No | None | N/A | Yes |
| M5 | Major | "first repeat to last" misdefines the hull | Apply | Applied | No | CLI ×4, `parameters_file/{mod,to_toml}.rs`, testdata, digests, specs | Pass | No |
| M6 | Major | "catalog never meets a contig end" | Apply | Applied | No | `region_typing/mod.rs`, specs | N/A (doc) | No |
| Mi1 | Minor | corrupt row panics in `BundleMember::new` | Defer | Deferred (doc updated) | No | `region_typing/mod.rs` | N/A | Yes |
| Mi2 | Minor | ASCII `expect` in `members_json` | Defer | Deferred | No | None | N/A | Yes |
| Mi3 | Minor | lookup `expect` without PANIC-FREE | Apply | Applied | No | `repeat_catalog/segments.rs` | Pass | No |
| Mi4 | Minor | `TheFiveFlags` | Apply | Applied | No | `cli/run_ground.rs` | Pass | No |
| Mi5 | Minor | right contig-end off-by-one | Apply | Applied | No | `region_typing/mod.rs` (test) | Pass | No |
| Mi6 | Minor | `i64::MAX` refusal untested | Apply | Applied | No | `psp/segmentation_section.rs` (test) | Pass | No |
| Mi7 | Minor | differential scan side `..default()` | Apply | Applied | No | differential test | Pass | No |
| Mi8 | Minor | dropped bundle not charged `FlankClamped` | Ask | Deferred | Yes | None | N/A | Yes |
| Mi9 | Minor | demotion of old parameters files undocumented | Apply | Applied | No | `str_bundle_loci.md`, `parameters_file.md` | N/A (doc) | No |
| Mi10 | Minor | garbled checksum note | Apply | Applied (with M5) | No | `cli/cross_platform_digests.rs` | Pass | No |
| Mi11–Mi17 | Minor | refactors and API shape | Defer | Deferred | No | None | N/A | Yes |

## 3. Questions asked and answers

1. **M2** — Drop a member wholly inside another at typing (recommended, after measuring), or let
   step 3's parse take overlapping members? — **Awaiting the owner.**
2. **Mi8** — Charge a bundle dropped at a contig end as `FlankClamped`, or state in §5.4 that the
   tally does not? — **Awaiting the owner.**

## 4. Per-finding log

### M1 — a windowed catalog read types a chain of short repeats as a bundle
- **Final status:** Applied.
- **Implementation:** `RepeatCatalogHeader::reach_into_window` takes the reader's bundle cap and
  pads by `max(longest tract, cap) + radius`; `rows_for` passes it. Doc of both says what the pad
  bounds (every member of a bundle at most the cap) and what it does not (a satellite's extent).
- **Test:** `a_window_over_one_end_of_a_bundle_types_the_whole_chain` — a 74-base, three-member
  chain whose longest member is 20 bases must come back from a read over its first member exactly
  as the scan types it; a 180-base, six-member chain must come back a satellite with the scan's
  start.
- **Verification:** the test passes; with `.max(max_str_bundle_len)` removed it fails (mutation run
  in the container, file restored and compared byte for byte). HG002's typing after the fix is
  identical to the implementer's step-2 typing (`diff` of every row: 0 lines).
- **Residual risk:** a satellite's extent from a windowed read can still be cut (out of scope
  observation in the review).

### M2 — members overlap and nest
- **Final status:** Applied with adaptation — the behaviour is the owner's decision; what was
  applied is the record. `RegionKind::SsrBundle`'s doc says the order is by detected span, the
  members are not disjoint, and gives the HG002 counts (551 nested, 734 overlapping, 520 where the
  last-listed member does not end the hull, 6 out of order). `str_bundle_loci.md` §5.1 gains a
  "found in review, open for step 3" paragraph naming what in §4, §6.1 and §7.1 it breaks and the
  three choices, with this review's recommendation; `typed_regions.md`'s type block no longer says
  "coordinate-ordered".

### M3 — the hull is not pinned to the detected spans
- **Final status:** Applied. `the_hull_ends_past_a_ragged_last_member`: `(TCTC)8`, 9 bases, `(CAG)8`
  then `CA`; asserts the member ends on 265 and the hull on 267.

### M4 — the catalog route's no-clean-cut fallback is untested
- **Final status:** Deferred. A fixture needs a repeat impure enough that `minimal_trim` fails yet
  clearing the prefilter and score floor, inside a bundle; finding one is trial and error against
  the detector, and the fallback itself was verified correct by reading. Carry forward.

### M5 — "first repeat to last"
- **Final status:** Applied. One wording everywhere — the hull runs from the first repeat's start
  to the furthest end of any of them: the three subcommands' help, `type-regions`' help (also
  dropping the backticked `Satellite`), `RepeatRouting::max_str_bundle_len`, the `to_toml.rs` note,
  `every_shape_as_written.toml`, the glossary, `SsrBundle`'s doc, and the specs' "385 … last
  member's" sentences (now "every member's"). The help also defines a bundle as "two or more
  repeats within 15 bases of each other" instead of "closer together than a clean flank".
- **Checksums:** the three parameters checksums moved, and **this rewording alone moves them**:
  restoring the earlier sentence in the file text each failing test printed gives the
  implementer's recorded `4e643855…`, `9265d027…`, `0cda09cc…` exactly. New values `a5b4997b…`,
  `4859eacd…`, `0573095b…`; the VCF checksums did not move. Notes rewritten (Mi10). The
  implementer's own claim — that reverting their three changes restores the original checksums —
  was not re-derived: the old census digests cannot be recomputed from the printed file.

### M6 — "the catalog never meets a contig end"
- **Final status:** Applied (doc): `resolve_features`' doc and both specs now say both routes apply
  the rule and why it rarely fires on the catalog's.

### Mi1 — a corrupt catalog row panics in `BundleMember::new`
- **Final status:** Deferred. The right fix validates spans in the catalog decoder, which predates
  this change and also underflows in `detected_interval`. `new`'s panic doc now says so.

### Mi2 — `members_json`'s ASCII `expect`
- **Final status:** Deferred: depends on whether the FASTA reader admits non-ASCII bases, not
  checked here; the same `expect` already exists for segments in the same file.

### Mi3, Mi4, Mi5, Mi6, Mi7, Mi9, Mi10
- **Mi3:** `PANIC-FREE` comment naming why every key is present and why a duplicate key is
  harmless. **Mi4:** `TheFiveFlags` → `TheCommandLineFlags`, its doc no longer counts flags.
  **Mi5:** a hull ending on base 199 of 200 (last base set to `C` so the detector cannot stretch
  onto it) stays a bundle. **Mi6:** `max_str_bundle_len_bp = u64::MAX` refused naming
  `segmentation.repeat-tract-criteria.max-str-bundle-len-bp`. **Mi7:** `scanned` now goes through
  `scan_config`, which names every field. **Mi9:** both specs say an old parameters file is demoted
  by a run whose census was built after. All validated by the runs in §1.

### Mi8, Mi11–Mi17
- **Deferred.** Mi8 is question 2. Mi11 (one flattened clap struct), Mi12 (one constant for the
  catalog's recorded cap), Mi13 (`Position` in `BundleMember`), Mi14 (moving `members_json`),
  Mi15 (grouping `resolve_features`' parameters) are refactors beyond a minimal fix and better taken
  with step 3's changes to the same types. Mi16 waits on M2. Mi17 is a one-line doc addition left
  with Mi12.

## 5. Deferred findings to carry forward
- M4, Mi1, Mi2, Mi8, Mi11–Mi17, as above.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no bench in `benches/` reaches `rows_for`; the one production change widens a read pad
only on contigs whose longest tract is shorter than the cap.

## 10. Commands run
All through `scripts/dev.sh` by absolute path; scripts under `tmp/review_2026-10-11_str-bundles-typing/`.
- `baseline.sh`: `cargo test --release --no-fail-fast --lib`; `--test ng_repeat_catalog_differential`
- `probe.sh`: the temporary M1 probe test (deleted after)
- `verify1.sh`, `verify2.sh`: touched-module filters (`region_typing repeat_catalog
  psp::segmentation_section segmentation_inputs parameters_file cli::typed_regions cli::run_ground
  cross_platform_digests`), the differential, and the M1 mutation
- `verify3.sh` (with `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks`): full library
  tests, `cargo fmt --check`, `cargo clippy`, the typing dump on HG002's regions and on
  `tmp/bundle_bench/runs/pcr_1s_30x` (outputs in `stutter_check/`; nothing under `runs/` written)
- `fmt.sh`: `cargo fmt` on three files, then the check

## 11. Command results
- baseline lib → 101, 4,995 passed, the 3 known release failures; differential → 0, 12 passed.
- probe → the window typed `bundle [200, 254]`, the whole contig `satellite [200, 380]`.
- verify2 → touched filters 0 (499 passed); differential 0 (13 passed); mutation → 101 (the new
  test fails), file restored.
- verify3 → lib 101 (4,996 passed, the 3 known); fmt clean after `fmt.sh`; clippy 101 on two
  untouched files; HG002 typing: 2,410 bundles, 328 satellites, 20,183 `ssr_locus`, 70,608
  `generic`, every row identical to `tmp/bundles_step2/typed.tsv`, and the `ssr_locus`/`generic`
  rows and all spans identical to main's `tmp/bundles/typed.tsv`; stuttering benchmark: 3,000 STR
  segments as planted, **500 of 500** bundles typed with the drawn members and motifs, 0 invented.

## 12. Notes
- `PROJECT_STATUS.md` was not updated: it has no block for the STR bundle work, and the branch's
  previous review (`str_bundles_simulator_2026-10-11.md`) did not create one either. The owner may
  want one before merge.
- Nothing committed, per the owner's instruction.
