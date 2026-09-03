# Fix Application Report: ng_psp_mode_b2_b3_2026-09-03.md

**Date:** 2026-09-03
**Source review:** `doc/devel/reports/reviews/ng_psp_mode_b2_b3_2026-09-03.md`
**Source state reviewed against:** 4fcdcd5b + the uncommitted B2/B3 diff, branch `ng-psp-mode`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 12
- Nits: 11 (grouped)

### Outcome totals
- Applied: 15 (M1–M3, Mi1–Mi12) plus 9 nits
- Deferred: 2 (see §5)
- Applied with adaptation / Already fixed / Disputed / Failed validation / Blocked / Superseded / Awaiting: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib 'ng::run::gatherer'` → 0, **16 passed** (same count: the new assertions
  landed inside existing tests rather than as new ones)
- `cargo test --lib 'ng::psp'` → 0, 417 passed; `cargo test --lib 'ng::read::input'` → 0, 228 passed
- The real-data oracle, re-run after the fixes with `NG_TWICE=1`: 183,807 records equal,
  header equal, second gather byte-identical (948,689 bytes), **exit 0** — identical figures
  to the pre-fix run, so nothing the fixes touched moved the measurement.
- Performance check → skipped: no `benches/`-reachable code changed.

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| M1 | Major | timestamp unasserted, clock-stamping gatherer passes | Applied | gatherer.rs |
| M2 | Major | env-count parse panics outside the exit contract | Applied | ng_psp_gather_oracle.rs |
| M3 | Major | failing dir entry silently drops a sample | Applied | ng_psp_gather_oracle.rs |
| Mi1 | Minor | empty-ground test cannot see a skipped segment | Applied | gatherer.rs |
| Mi2 | Minor | example never checks the tract count it prints | Applied | ng_psp_gather_oracle.rs |
| Mi3 | Minor | homopolymer fixture blind to symmetric defects | Applied (doc) | gatherer.rs |
| Mi4 | Minor | chain-equality pin covers only generic ground | Applied | gatherer.rs |
| Mi5 | Minor | `NG_TWICE` absent from usage | Applied | ng_psp_gather_oracle.rs |
| Mi6 | Minor | a wasted third `open()` | Applied | ng_psp_gather_oracle.rs |
| Mi7 | Minor | `from_ref` of a temporary | Applied | ng_psp_gather_oracle.rs |
| Mi8 | Minor | zero counts blame the wrong input | Applied | ng_psp_gather_oracle.rs |
| Mi9 | Minor | segmentation scaffold at three copies | Applied | gatherer.rs |
| Mi10 | Minor | three shapes of "open at the unusual settings" | Applied | gatherer.rs |
| Mi11 | Minor | `gather_and_compare` at ~106 lines | Applied | ng_psp_gather_oracle.rs |
| Mi12 | Minor | four naming repairs incl. the magic zstd key | Applied | ng_psp_gather_oracle.rs, gatherer.rs, psp/writer.rs, psp/mod.rs |
| Nits | Nit | 9 of 11 applied; 2 deferred | Applied / Deferred | both files |

## 3. Questions asked and answers
None — non-interactive. The review's two open questions were resolved as recorded there:
a tract-free walk is now a refusal (Q1), and the `u32` contig-length narrowing is left to the
upstream `regions` layer (Q2).

## 4. Per-finding log (compressed)

- **M1 — Applied.** `gathering_the_same_sample_twice_gives_identical_bytes` now reopens the
  first file and asserts `header.writer.created == provenance().created`, with a comment
  saying why the byte comparison alone cannot see a self-stamping clock. This is the probe
  that killed the reviewer's surviving mutation A.
- **M2 + Mi8 — Applied.** `how_many` returns `Result`, names the variable and the offending
  value, and refuses zero there; both call sites propagate with `?`.
- **M3 — Applied.** The directory listing collects fallibly (`Result<Vec<_>, _>`), so an
  entry error ends the run instead of removing a sample from the oracle.
- **Mi1 — Applied.** The empty-ground test asserts `regions_handled == 1` beside
  `regions_in`, with the message spelling out that the psp is byte-identical whether the
  segment was walked or skipped, so only the tally can tell.
- **Mi2 — Applied.** A walk with zero tract records now prints `!!` and returns `Ok(false)`,
  naming what to widen.
- **Mi3 — Applied.** The tract test's doc records the homopolymer limitation, the measured
  fact behind it (a flank-swap mutant passed), and where the class *is* pinned.
- **Mi4 — Applied.** `direct_walk` takes a segmentation, and the tract test now compares its
  walk against the bare direct-mode chain over tract ground — closing the gap in the warrant
  the example's doc cites.
- **Mi5 — Applied.** `NG_TWICE=1` is in the usage string.
- **Mi6 — Applied.** The header is cloned from the first gatherer before `write_psp` consumes
  it; three opens become two, and the clone is no longer billed to the comparison walk.
- **Mi7 — Applied.** `let alignments = [path.to_path_buf()];` hoisted above the closure.
- **Mi9 + Mi10 — Applied.** `build_segmentation(segments, analysed)`, `fixture_bounds()`,
  `generic_segment(contig, length)` and `open_gatherer_over(...)` land in the test module;
  the three scaffold copies and the three open shapes each collapse to one.
- **Mi11 — Applied.** `header_matches`, `gathered_twice_is_byte_identical` and
  `file_matches_walk` are named functions the module doc's "what each comparison proves" list
  now maps onto; `gather_and_compare` keeps the gather, the timing and the summary.
- **Mi12 — Applied.** `header_the_gatherer_fixed`, the closure `open_gatherer`,
  `chr2_only_ground`, and a new exported `ZSTD_COMPRESSION_LEVEL_KEY` in `psp::writer` that
  the writer's own two sites and the example all use.
- **Nits — Applied:** `NG_TWICE` now requires the value `1`; `report_first_difference`
  exhaustively destructures the record so a new field cannot be misreported; a header
  mismatch names the differing field group; the doc's "exits non-zero on the first sample"
  becomes "checks every sample"; `work` → `work_dir`; `read_back` → `records_read_back`;
  `first`/`second` → `first_psp`/`second_psp`; "re-walk" → "comparison walk";
  `first_regions_of`'s shadowing parameter → `count`; the `expect` carries a `PANIC-FREE:`
  note. **Deferred:** the two `mut`-bindings-during-init nits (the shape matches the
  pre-existing tests beside them, so changing one of three is churn).

## 5. Deferred findings to carry forward
- **`entry.length as u32`** — the `ContigBounds`/`regions` layer is `u32` repo-wide and the
  sibling example does the same; an upstream question, not this diff's.
- **The run-level shared test-fixture module** — still deferred from B1; Mi9's local helpers
  are the interim it will absorb.

## 6–8. Disputed / Failed validation / Blocked
None.

## 9. Performance check
Skipped — no `Apply` touched code reachable from `benches/`.

## 10–11. Commands run and results
In the container via `scripts/dev.sh`, on the tree committed: `cargo fmt` / `--check` → 0;
`cargo clippy --all-targets --all-features -- -D warnings` → 0; `cargo test --lib
'ng::run::gatherer'` → 0 (16); `'ng::psp'` → 0 (417); `'ng::read::input'` → 0 (228);
`NG_TWICE=1 cargo run --release --example ng_psp_gather_oracle -- …` → exit 0 with the
figures quoted in §1.

## 12. Notes
The zstd-key constant is the one fix that reaches outside the reviewed diff (two call sites
in `psp/writer.rs` plus a re-export). It changes no behaviour — the literal it replaces is
byte-identical — and `ng::psp`'s 417 tests cover the writer's use of it.
