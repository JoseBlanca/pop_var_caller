# Fix Application Report: str_bundles_calling_2026-10-11.md

**Date:** 2026-10-11 (runs finished 2026-10-12)
**Source review:** `doc/devel/reports/reviews/str_bundles_calling_2026-10-11.md`
**Source state reviewed against:** `str-bundles` at 8fee2920 plus the uncommitted step-4 diff
**Execution mode:** non-interactive (session lead's brief)
**Overall status:** Partial — the Blocker and every correctness Major fixed; test-only and refactor findings deferred

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 12
- Minors: 19
- Nits: grouped

### Outcome totals
- Applied: 6 (B1, M1, M2, M5, M6, M7)
- Applied with adaptation: 0
- Deferred: 6 Majors (M3 awaiting the owner, M4, M8–M12), 17 Minors
- Applied Minors: Mi1 (and the stale texts listed under M6)
- Disputed / failed / blocked / superseded: 0

### Validation summary
- `cargo fmt -- --check` → clean.
- `cargo clippy --release --all-targets` → 2 warnings, both pre-existing in untouched code.
- `cargo test --release --no-fail-fast --lib` → 5,064 passed, 3 failed (the known debug-only
  `should_panic` tests), 11 ignored.
- `cargo test --release --no-fail-fast -j 2 --tests --bins` → integration targets all pass.
- Cross-platform VCF checksums (`cli::cross_platform_digests` tests, in the lib run) → pass,
  unchanged: their fixture holds no bundle.
- `cargo doc`, `cargo audit`, benches → not run (no dependency change; no bench covers bundles).

### Unresolved high-priority findings
- M3 — the part-repeat regime at bundles whose periods share no divisor above 1: owner decision.
- M4 — no test from selection to a VCF record at a bundle (the benchmark exercises it).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | Edge deletion written twice / its reads left out | Apply | Applied | No | `src/locus_generation/ssr.rs` | Tests + 4 benchmark runs | Residual risk (HG002, 2 regions) |
| M1 | Major | Zero-copy run costs 0.88 a read | Apply | Applied | No | `src/calling/likelihood/ssr_emission.rs` | Tests | No |
| M2 | Major | Part-repeat place in a run too short | Apply | Applied | No | same | Tests | No |
| M3 | Major | No part-repeat regime at gcd-1 bundles | Ask | Deferred | Owner | None | N/A | Yes |
| M4 | Major | No bundle selection-to-VCF test | Defer | Deferred | No | None | N/A | Yes |
| M5 | Major | Misleading 1.x psp refusal | Apply | Applied | No | `src/psp/mod.rs`, `src/psp/header.rs` | Tests | No |
| M6 | Major | Stale help text and docs | Apply | Applied | No | `src/cli/command_line.rs`, `src/cli/call_from_alignments.rs`, `src/calling/likelihood/ssr_emission.rs`, `src/calling/mod.rs`, `src/calling/inference/repeat_tract_parameters.rs`, `doc/devel/ng/spec/locus_generation_ssr.md` | Tests (renamed panic text) | No |
| M7 | Major | Partial exclusion relies on row order | Apply | Applied | No | `src/run/callers.rs` | Tests + benchmark | No |
| M8 | Major | `gather_for_bundle` untested at ≥2 runs | Defer | Deferred | No | None | N/A | Yes |
| M9 | Major | `fill_emissions` flank shift untested | Defer | Deferred | No | None | N/A | Yes |
| M10 | Major | No bundle in the digest fixture | Defer | Deferred | No | None | N/A | Yes |
| M11 | Major | Stand-in segment fields at bundles | Defer | Deferred | No | None | N/A | Yes |
| M12 | Major | Non-discriminating selection/prior tests | Defer | Deferred | No | None | N/A | Yes |
| Mi1 | Minor | All sets pruned → part-repeat | Apply | Applied | No | `src/calling/likelihood/ssr_emission.rs` | Tests | No |
| Mi2–Mi19 | Minor | see review | Defer | Deferred | No | None | N/A | Yes |

## 3. Questions asked and answers

1. **M3** — keep §7.3's rule (part-repeat only when no slip set reaches the change) and document
   how many bundles it covers, or add each non-whole-unit run's part-repeat place to the slip-set
   sum? **Answer:** awaiting the owner.

## 4. Per-finding log

### B1 — A deletion across a hull's outer edge is written twice, or its reads are left out
- **Initial decision / final status:** Apply / Applied.
- **Why bundles missed it** (the lead's questions): the own-alignment route *is* used at bundles
  (same `classify_read`), and gathering *does* see bundle records. What fails is upstream of both:
  in the three bundles concerned the mapper never writes a deletion crossing the hull's edge — it
  writes substitutions (`bun00119`, `150M`) or a flank deletion plus an insertion at the hull's
  first base (`bun00018` `64M4D2M6I78M`, `bun00031` `68M1D2M2D80M`) — so no read carries a flank
  deletion, and the flat-gap delimiter either spells the flank change into the hull (written
  twice, once by each path) or cannot anchor the flank (every carrier read a partial, which
  calling leaves out at a bundle).
- **Implementation:** `MeasuredRepeat` gains `ground: MeasuredGround { Segment, BundleHull }`;
  in `classify_read`, after the delimiter and only at a `BundleHull`,
  `flank_deletion_the_delimiter_implies` turns a read whose flank the delimiter shortened (or
  could not anchor, with the CIGAR deleting bases in that flank) into a complete observation:
  the flank's longest stretch matching from its far end (≥ 5 bases), the deleted bases, the
  spelling, in smallest form. No 8-base overhang rule; no parameter.
- **Tests added:** `a_bundle_read_whose_delimiter_deletes_flank_bases_records_them_as_a_flank_deletion`
  (the `bun00119` window, and the reference read unchanged),
  `a_bundle_read_with_a_deletion_inside_the_flank_records_it_as_a_flank_deletion` (the `bun00018`
  window, `10M4D2M6I61M`). Both use the benchmark's own sequences.
- **Probe:** on new psps of three samples the carrier reads became complete observations with
  `(3,0)`, `(5,0)` and `(5,0)` flank deletions and exactly the truth's spellings.
- **Measured** (new psps under `tmp/review_fix/bench/`, true stutter; bundle cells with a variant
  right / wrong / missed, false calls):

  | run | before review | after review |
  |---|---|---|
  | 1 sample, PCR, 30× | 205 / 3 / 0, 0 false | 207 / 1 / 0, 0 false |
  | 1 sample, PCR, 5× | 98 / 30 / 80, 0 false | 98 / 30 / 80, 0 false |
  | 50 samples, PCR, 30× | 10,349 / 51 / 35, 3 false | 10,374 / 26 / 35, 3 false |
  | HG002, 30× (per region) | 824 / 158 / 149, 32 false | 822 / 160 / 149, 32 false |

  At 50 samples the 374 bundle cells with an edge deletion went from 347 right, 24 wrong and 3
  "lengths right" to 372 right and 2 wrong. Segment cells: 3,000, 3,000 and 150,000, none differ
  from main.
- **Residual risk (HG002, 2 of 1,131 regions right → wrong):** chr4:28,575,007 — the record now
  spells the truth's flank SNP `G>C` inside a widened bundle allele; the sequence is the truth's,
  and the region scorer does not credit a record widened past the region. chr3:81,721,205 — the
  haplotype carrying a one-base deletion just past the hull has two crossing reads; one now
  carries the flank deletion and the other does not, so each key has one read and the allele
  misses the two-read bar. Reported to the owner (review §4, question 2).

### M1 — A run with no copy costs 0.88 on every read
- Applied. `bundle_emission`: probability `[1.0]` for a zero-copy run; `part_repeat_emission`
  skips zero-copy runs in the product of the other runs' shares. Tests:
  `a_run_with_no_copy_contributes_a_factor_of_one` (new, exact value);
  `a_run_with_no_copy_in_the_allele_cannot_slip` now asserts the exact part-repeat value instead
  of an upper bound. Test-first not run separately: the old code returns 0.88 × the expected
  value, which the new assertion's 1e-12 relative tolerance rejects. No benchmark cell moved.

### M2 — Part-repeat place in a run too short for the change
- Applied. A place needs a copy and `stretch + Δ ≥ 0`. Test
  `a_run_with_no_copy_is_not_a_place_for_a_part_repeat_change` (would panic in debug before, and
  return half the value in release). The one-run regression still passes bit for bit.

### M3 — gcd-1 bundles have no part-repeat regime
- Deferred: an owner decision (question 1).

### M4 — No bundle test from selection to a VCF record
- Deferred: it needs a cohort observation fixture through `call_one_ssr_bundle_locus`,
  `evidence_for_output`, `assemble_record` and `encode`. Checked instead: the run counts are read
  off the post-drop table, and in the 50-sample VCF all 452 bundle records carry `REPCNB` with one
  value per genotype allele and no `.` for a called allele.

### M5 — Misleading 1.x psp refusal
- Applied. `PspReadError::OlderVersion { path, found, oldest_major }`:
  "x.psp is format 1.2, written by an older pvc; this reader reads format 2.0 and later, so write
  the psp again with this pvc". `UnsupportedVersion` keeps newer majors. Tests:
  `a_format_1_header_and_a_major_before_it_are_refused_as_versions` now expects `OlderVersion`;
  the message test in `psp/mod.rs` covers the new text. Public error enum gained a variant — no
  exhaustive match outside `psp` existed (the crate compiled unchanged elsewhere).

### M6 — Stale help text and docs
- Applied: `call-from-alignments --help`, the `call_from_alignments.rs` module doc, the
  `censored_emission` comment, `calling/mod.rs`'s doc, `tract_candidates`' doc, panic message and
  its `should_panic` test (now expects "this one-tract table does not score a bundle"), and
  `locus_generation_ssr.md`.

### M7 — Partial exclusion relies on row order
- Applied: a release `assert!` in `call_one_ssr_bundle_locus` that no complete row follows the
  complete prefix. No benchmark run tripped it.

### M8–M12, Mi2–Mi19 — Deferred
- Test-only findings (M8, M9, M10, M12, Mi3, Mi4, Mi17) and refactors (M11, Mi7, Mi10–Mi15) left
  for a test-and-cleanup pass; none changes a number the benchmarks show. Mi5, Mi6, Mi8 and Mi9
  are behaviour choices for the owner; Mi2 needs either the claim narrowed or the slip floor taken
  from the run's segmentation; Mi18 and Mi19 are noted in the spec's cost line and the review.

### Mi1 — All slip sets pruned
- Applied: `any_slip_set_reaches` decides whether Δ was reached; if it was and every set was
  pruned, the read scores 0 rather than the part-repeat value.

## 5. Deferred findings to carry forward
- M3 (owner), M4, M8, M9, M10, M11, M12; Mi2–Mi19.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check

- **Triggered:** no criterion bench covers the bundle path, so measured directly instead.
- **Cost of a bundle locus** (one thread, 50 samples, PCR, 30×, psps over the 500 bundle hulls
  ±60 bases): calling takes 1.98 s with this build against 0.53 s with step 3's, which sets
  bundles aside — **1.45 s for 500 bundles, about 2.9 ms each**. Over 2,000 segment regions the
  calling pass takes 4.0 s, generic loci in their flanks included, so **at most 2 ms a segment**.
  On the whole 50-sample run (8 threads) the calling pass is 6.8–8.1 s across builds; bundles
  are about 0.2 s of it. The perf sub-agent's candidates (per-slip-set `ln` calls, per-set
  allocations, per-read recomputation of run segmentation and slip probabilities) are recorded in
  `tmp/review_2026-10-11_str-bundles-calling/float_perf.md` with a measurement plan; none applied.

## 10. Commands run
- `scripts/dev.sh bash tmp/review_fix/build.sh` (`cargo build --release --bin pvc`)
- `scripts/dev.sh bash tmp/review_fix/test.sh <filter>` (`cargo test --release --lib <filter>`)
- `scripts/dev.sh bash tmp/review_fix/test_all.sh` (lib; `--tests --bins -j 2`)
- `scripts/dev.sh bash tmp/review_fix/lint.sh` (fmt check, clippy), `tmp/review_fix/fmt.sh`
- `scripts/dev.sh bash tmp/review_fix/bench/run_all.sh` (edge fix, then final binary as `x_final`)
- `DEV_EXTRA_MOUNT=… scripts/dev.sh bash tmp/step4/hg002_two_binaries_30x.sh tmp/review_fix/hg002 tmp/bins/pvc_step4_prereview tmp/bins/pvc_step4_reviewed`
- `scripts/dev.sh bash tmp/review_fix/cost.sh`; probes `tmp/review_fix/probe*.sh`, `hg2probe.sh`
- `scripts/dev.sh bash tmp/review_fix/cleanup.sh` (all generated psps deleted)

## 11. Command results
- lib tests → exit 101, 5,064 passed / 3 known failures.
- integration tests → all pass.
- fmt → clean; clippy → 2 pre-existing warnings.
- benchmarks → table in B1.

## 12. Notes
- The final binary is `tmp/bins/pvc_step4_reviewed`; the pre-review build is
  `tmp/bins/pvc_step4_prereview`.
- `PROJECT_STATUS.md`: a block for STR bundles was added pointing at these two reports; the
  "Current focus" lines were left to the session lead.
- The spec records all of this: `str_bundle_loci.md` §6.1 (the new measurement rule), §7.3 (the
  three emission corrections), §10 step 4 (after-review table, the two HG002 flips, cost).
