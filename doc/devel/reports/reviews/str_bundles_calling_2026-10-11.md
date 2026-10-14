# Code Review: str_bundles_calling
**Date:** 2026-10-11
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step 4 of `str_bundle_loci.md` — STR bundles called as one STR locus (uncommitted diff on `str-bundles`, HEAD 8fee2920)
**Status:** Request-changes (the Blocker and the three correctness Majors are fixed in the same session; see [fixes_applied_str_bundles_calling_2026-10-11.md](fixes_applied_str_bundles_calling_2026-10-11.md))

---

## 1. Scope

- **What was reviewed:** the uncommitted diff of branch `str-bundles` over 8fee2920 (31 files,
  +2,890 / −292) and the untracked `src/region_typing/run_parse.rs`. Design:
  [`str_bundle_loci.md`](../../ng/spec/str_bundle_loci.md) §7 and its "Changed during
  implementation" notes; implementation report
  [`str_bundles_calling_2026-10-11.md`](../implementations/str_bundles_calling_2026-10-11.md).
- **Reviewed against:** 8fee2920 plus the working tree as found.
- **In scope:** run parse; `select_ssr_bundle` and the shared selection passes; `bundle_emission`,
  slip sets, the part-repeat regime; `TractScoringFits::gather_for_bundle`, `bundle_candidates`;
  `fill_ssr_bundle_seed`; `call_one_ssr_bundle_locus` and `LocusEvidence::Ssr { shape }`; VCF
  (`RU`/`PERIOD` `Number=.`, `REPCNB`, no-run bundles); psp `OLDEST_READABLE_MAJOR = 2`; partials
  at bundles; the re-recorded cross-platform checksums. Also the bundle read measurement in
  `locus_generation/ssr.rs`, because the edge-deletion defect lives there.
- **Out of scope:** the periodicity gate at bundles (the owner's pending decision); the typing and
  psp payload of steps 2–3 (reviewed earlier).
- **Categories dispatched** (read-only sub-agents in this worktree, no worktrees and no builds,
  per the owner's rule on disk; so no mutation was run by a sub-agent — every "would survive" below
  is reasoned):
  - reliability, split in two (the parse and emission core; selection, prior, driver and VCF) —
    always;
  - errors + defaults + extras — always, and "diff matches stated intent";
  - naming + idiomatic + smells + refactor safety + module structure — always / multi-file;
  - float portability + a brief performance pass — float arithmetic reaches every genotype.
  Per-category files: `tmp/review_2026-10-11_str-bundles-calling/`.
- The lead investigated the edge-deletion defect directly (B1) with read-level probes.

## 2. Verdict

Request-changes. One Blocker (the edge-deletion double description, 24 of the 51 wrong bundle
cells at 50 samples) and two Majors that put wrong numbers into the bundle likelihood. None of the
three moves a segment.

## 3. Execution status

Run by the lead, in the dev container, after the fixes (pre-fix numbers come from the
implementer's runs and were re-derived from their scored files):

- `cargo test --release --no-fail-fast --lib` → exit 101: **5,064 passed, 3 failed, 11 ignored**;
  the three failures are the known debug-only `should_panic` tests
  (`an_offset_past_the_reference_panics` ×2, `lgamma_panics_on_non_positive_in_debug`).
- `cargo test --release --no-fail-fast -j 2 --tests --bins` → the 11 integration targets pass;
  the only failures are the same three lib tests.
- `cargo fmt -- --check` → clean after formatting the touched files.
- `cargo clippy --release --all-targets` → two warnings, both pre-existing in code this work did
  not touch (`locus_generation/ssr.rs` tally's `type_complexity`; `run/edge_deletion_gathering.rs`
  `single_match`).
- Not run: `cargo doc`, `cargo audit` (no dependency change), `cargo bench` (no bench covers the
  bundle path).
- Findings labelled "Needs verification": 0. Every "a mutation would survive" claim is reasoned,
  not run.

## 4. Open questions and assumptions

1. **(M3)** At bundles whose runs' periods have no common divisor above 1 — any homopolymer run, or
   periods 2 and 3 — every length change is reachable by some slip set, so the part-repeat regime
   never runs and a one-base sequencing indel is priced only as two runs slipping in opposite
   directions (≈ 0.0023 against a segment's ≈ 0.0095 under the shipped constants). The code
   follows §7.3. Should the part-repeat place be added to the slip-set sum for runs where the
   change is not whole units, or is the spec's rule kept and documented with a count of the
   bundles it covers?
2. **(Residual risk of B1's fix)** The new rule records a flank deletion only where the delimiter
   itself shortened or failed to anchor a flank. On HG002 it moved two regions from right to
   wrong (§10 of the spec, step 4 after review). Is a one-read-each split of one haplotype into
   two keys acceptable at that rate, against 25 cells gained at 50 samples?

## 5. Top 3 priorities

1. **B1** — a deletion across a hull's edge written twice, or its reads all dropped: the largest
   source of wrong bundle cells. Fixed.
2. **M1** — a run the allele holds no copy of scales every read by 0.88. Fixed.
3. **M4** — no test takes a bundle from selection to a VCF record; four wiring mutations survive.
   Deferred (the benchmark exercises the path; a driver fixture is the next task).

## 6. Findings

### Blocker

#### B1: [src/locus_generation/ssr.rs](../../../../src/locus_generation/ssr.rs#L1236) — A deletion across a hull's outer edge is written twice, or its reads are all left out
- **Confidence:** High (probed read by read on three bundles; measured).
- **Problem:** Bundles take the same `classify_read` as segments, so the own-alignment route
  (a CIGAR deletion crossing the edge → flank deletion in smallest form) is used at bundles. It
  rarely fires there: the simulator lays an edge deletion on a length allele, so the read is often
  the reference's length or longer, and the mapper writes the event as substitutions or as a
  deletion inside the flank plus an insertion at the hull's first base. Concretely:
  - `bun00119` (`CTCCGA|(AAAC)6`, truth `CTC|AAC(AAAC)6`, same length): reads are `150M`. The
    flat-gap delimiter, which prices nothing inside the hull and treats a gap beside the junction
    softly, deleted two flank bases and spelled `AC` + the reference hull (50 bases). The generic
    path called `CCGA→CAAC` beside it. Two records, one event.
  - `bun00018` (`…TCCGA|(AAAT)6`, truth `…|AAT(AAAT)7`): reads `64M4D2M6I78M`. No deletion
    crosses the edge; the delimiter cannot anchor a flank whose last five bases are gone, so
    every carrier read was a one-sided partial (11 reads covering the whole hull in `sim009`),
    which calling leaves out at a bundle. The generic path called the flank deletion alone.
  - `bun00031`: the same as `bun00018`, written `68M1D2M2D80M`.
  Gathering does see bundle records; it had nothing to gather, because no bundle allele carried a
  flank deletion. Segments escape because their delimiter prices an off-unit insertion inside
  the repeat, so it keeps the flank and leaves the change to the generic path once; and a segment
  uses its partials.
- **Effect:** 24 wrong cells and 3 "lengths right only" among the 374 bundle cells with an edge
  deletion (50 samples, PCR, 30×).
- **Fix (applied):** at a bundle only, where the delimiter itself shortened a complete read's
  flank or could not anchor a side whose CIGAR deletes bases in that flank, read the window as the
  flank's longest stretch matching from its far end (at least five bases), the deleted bases, and
  the spelling; put it in smallest form. Guards: the read's own alignment covers the window and
  reaches five aligned bases beyond each edge; a side the delimiter anchored in full keeps its
  junction.

### Major

#### M1: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L1041) — A run the allele holds no copy of costs its model's no-slip share on every read
- **Categories:** errors (filed Blocker), reliability.
- **Confidence:** High (mechanism); measured effect: no benchmark cell moved.
- **Problem:** a zero-copy run gets range `(0, 0)` and probability `model.probability(0, period)`,
  and `gather_for_bundle` gives such a run the shipped model, whose no-slip share is 0.88. Every
  slip set, and every part-repeat place (through the other runs' `same_length_share`), is scaled
  by it. A sample with N reads of an allele that lost a whole run loses 0.88^N against genotypes
  without it: 0.02 at 30 reads. `for_bundle`'s unreachable mass counts nothing lost.
- **Fix (applied):** probability `[1.0]` for a zero-copy run; skip it in the part-repeat product.
  Test `a_run_with_no_copy_contributes_a_factor_of_one`; the old upper-bound test now asserts the
  exact value.

#### M2: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L1233) — The part-repeat regime places an indel in a run too short to take it
- **Categories:** errors, reliability.
- **Confidence:** High.
- **Problem:** a zero-copy run, or a stretch shorter than the contraction, counted as a place;
  `resize_at_the_end` saturates, so the rendered allele has another length than the read.
  Debug build: `substitution_probability`'s `debug_assert_eq!` panics. Release: that place scores
  0 and still divides the sum. Input: hull `GG (CA)8 GT (AAAG)4 C`, allele without the `CA` run,
  read one base short in `AAAG` → half the correct value.
- **Fix (applied):** a place needs a copy and `stretch + Δ ≥ 0` (never excluded at a segment, whose
  stretch is the whole allele). Test `a_run_with_no_copy_is_not_a_place_for_a_part_repeat_change`.

#### M3: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L1041) — No part-repeat regime at bundles whose periods share no divisor above 1
- **Confidence:** Medium (follows the spec; the effect is reasoned from shipped constants).
- **Problem and fix:** see open question 1. Owner decision.

#### M4: [src/run/callers.rs](../../../../src/run/callers.rs#L1571) — No test takes a bundle from selection to a VCF record
- **Categories:** reliability.
- **Confidence:** High that the tests are missing.
- **Problem:** nothing calls `call_one_ssr_bundle_locus`, the `SsrBundle` arms of
  `evidence_for_output`/`filter_for` (`run/records.rs`), or the `REPCNB` parse in
  `assemble_record`. Four wiring mutations would survive: parsing run counts from the table before
  the drop, returning no annotation for a bundle (gathering would then treat it as a generic
  record), putting bundles back in `filter_for`'s `Pass` arm, removing the partial slice.
- **Checked by the lead instead:** the code reads the post-drop table (`vcf/assemble.rs`), and on
  the 50-sample run every one of 452 bundle records has `REPCNB` with one entry per allele of the
  genotype and none `.` for a called allele.
- **Fix:** the driver-level test the reliability file specifies (two-run bundle, a candidate no
  genotype names, three partials; assert FORMAT, `REPCNB` after the drop, `STR;RU;PERIOD`, `DP`
  with partials and `AD` without; a second case for `tooManyAlleles`).

#### M5: [src/psp/header.rs](../../../../src/psp/header.rs) — A 1.x psp is refused with "this reader understands up to 2.0"
- **Categories:** errors.
- **Problem:** the message reads as though the file should load, and the docs promised "written by
  an older pvc; write it again".
- **Fix (applied):** a new `PspReadError::OlderVersion` naming the remedy; tests updated.

#### M6: [src/cli/command_line.rs](../../../../src/cli/command_line.rs) and [src/cli/call_from_alignments.rs](../../../../src/cli/call_from_alignments.rs) — `--help` and the module doc say bundles are set aside
- **Categories:** naming/smells, errors.
- **Fix (applied):** both describe bundle calling and `REPCNB`. Also fixed the stale
  `censored_emission` comment, `locus_generation_ssr.md`'s `ended_inside_a_bundle`, and the "nothing
  scores a bundle" doc, panic message and test in `repeat_tract_parameters.rs` / `calling/mod.rs`.

#### M7: [src/run/callers.rs](../../../../src/run/callers.rs#L1571) — Leaving partials out of calling depends on an unstated row order
- **Categories:** naming/smells, errors (Minor there).
- **Problem:** `partition_point` assumes `rebuild` writes complete rows first; if that changed,
  partials would reach `censored_emission` (0 in release, silently).
- **Fix (applied):** a release `assert!` that nothing complete follows the prefix.

#### M8: `TractScoringFits::gather_for_bundle` and `bundle_candidates` are untested with two or more runs
- **Categories:** reliability. Mean vs sum/first/max of the runs' substitution rates, a transposed
  per-run model table and per-run fallback counting all pass today. **Fix:** the test the
  reliability file specifies. Deferred.

#### M9: `fill_emissions`' left-flank shift at bundles has no test through `fill_emissions`
- **Categories:** reliability. Now load-bearing: B1's fix creates bundle alleles with flank
  deletions. The arithmetic was checked by reading and is right. **Fix:** the test with an
  interrupted first run. Deferred.

#### M10: The cross-platform and thread-count checks have no bundle
- **Categories:** reliability, float portability (Minor there). The fixture's checksums moved only
  with the header. **Fix:** add a two- and a three-run bundle to the fixture. Deferred.

#### M11: Segment-only fields hold stand-in values at a bundle
- **Categories:** naming/smells, errors. `SsrCandidate { bundle: Some(..), repeat_count = 1 }`, a
  context motif of `"A"` or the first run's, `TractAnnotation::motif()`/`period()` returning the
  first run's values. A second `SsrEmissionModel` would score bundles as one-repeat segments.
  **Fix:** an enum for the segment/bundle choice. Deferred (refactor across the calling seam).

#### M12: Two selection and prior tests cannot discriminate
- **Categories:** reliability. The selection fixtures never trigger the ±1 rescue on vectors and
  all alleles differ in length (`(11,3)` against `(8,5)`, both 31 bases, would discriminate);
  `fill_bundle_seed_share_per_candidate` is tested only on a flat prior, where deleting the
  normalisation passes. Deferred.

### Minor

- **Mi1** `collect_slip_sets` — when pruning removes every set reaching Δ, the read fell to the
  part-repeat regime. Fixed (`any_slip_set_reaches`).
- **Mi2** `run_parse.rs` tie to the reference count — `CACACACA` against hull `CACACAC` parses to 3
  copies, so the one-run bundle and the segment differ on the deepest contraction and on the
  stutter stratum. The "bit for bit" claim is narrower than stated. Deferred (narrow the claim or
  take the floor from the run's own segmentation).
- **Mi3** `run_parse.rs` copy cap unproven under the tie rule; no brute-force property test.
- **Mi4** pruning from three runs on has no test.
- **Mi5** a candidate with no copy in any run takes the default substitution rate and `Defaulted`
  at a bundle that has runs; the spec covers only bundles with no run.
- **Mi6** a flank length change at a bundle is priced with the first run's model, even at the right
  edge, and not counted as defaulted.
- **Mi7** `gather_for_bundle` zips each parse's runs with the bundle's without checking the counts.
- **Mi8** user-set `--discovery` and slippage re-fit are ignored at bundles without a word in the
  report.
- **Mi9** `RU`/`PERIOD` became `Number=.` for every VCF, segments included (spec §7.5 chose this;
  htslib readers now return a list).
- **Mi10** a missing `REPCNB` entry is written `.` silently; `with_run_counts_of_alleles` uses a
  struct-update spread and `run_counts_of_allele` a wildcard.
- **Mi11** `call_one_ssr_bundle_locus` copies about 35 of `call_one_ssr_locus`'s 50 lines;
  `gather_for_bundle` repeats a 15-field reset; the seed's spectrum lookup is copied — each a place
  the one-run equality can drift.
- **Mi12** `format_keys` and `one_sample_column` decide separately which of three cases a record
  is in.
- **Mi13** `run_parse` lives in `region_typing`, which never uses it.
- **Mi14** the ladder's keys are still called repeat counts, and "ties to the shorter" means "the
  vector that sorts first" at a bundle.
- **Mi15** two `let … else { unreachable!() }` on `SsrLocusShape` in `summarise_condition.rs`.
- **Mi16** `is_repeat_tract`'s doc says `STR` is always written when true; not for a no-run bundle.
- **Mi17** tests: the grid test's periods 2 and 3 cannot tell the gcd from "step 1"; no header test
  lists `REPCNB`; no test that bundles skip the slippage round; the bundle depth test has no
  partial.
- **Mi18** spec §8 says "a few terms per pair of runs"; two equal-period runs walk 21 sets, three
  331 before pruning.
- **Mi19** robustness: the parse allocates `states × allele_len` cells (16 bytes each), quadratic
  in allele length — gigabytes for a 10 kb long-read insertion. Irrelevant at 150-base reads.

### Nits

`build_template`'s doc names a parameter it does not take; parallel `run_ends`/`run_follows_a_run`
vectors; `TractOutcomes::bundles` naming; `let mut explanation = explanation`; `path_word` and the
count assertion say "repeat counts" for run parses; `set_count`'s redundant no-run branch;
`NEGLIGIBLE_SLIP_SET_PROBABILITY` could be `pub(crate)`; `&observation.alleles[0]` indexed before
the descriptive assert; the bare `runs.len() >= 3`; `unreachable_mass` not bit-identical to the
segment's for one run (reported only).

## 7. Out of scope observations

- `run/records.rs` `filter_for` ends in `_ => FilterVerdict::Pass` (pre-existing); bundles now pass
  through it, so a new verdict would be written as PASS. Separate change.

## 8. Missing tests to add now

Grouped by function; full specifications in the per-category files.

- `call_one_ssr_bundle_locus`: `a_two_run_bundle_is_called_and_written_with_run_counts_after_the_drop` (M4).
- `gather_for_bundle`: `gather_for_bundle_takes_the_mean_rate_and_one_model_per_run_per_cell` (M8).
- `fill_emissions`: `a_bundle_read_with_a_wider_left_flank_deletion_shifts_runs_by_the_flank_bases_given` (M9).
- `parse_into_runs`: `parse_into_runs_matches_brute_force_over_count_vectors` (Mi3).
- `bundle_emission`: `pruning_at_three_runs_changes_the_emission_by_less_than_the_threshold_times_the_sets` (Mi4);
  added in this session: `a_run_with_no_copy_contributes_a_factor_of_one`,
  `a_run_with_no_copy_is_not_a_place_for_a_part_repeat_change`.
- `classify_read` at a bundle — added in this session:
  `a_bundle_read_whose_delimiter_deletes_flank_bases_records_them_as_a_flank_deletion`,
  `a_bundle_read_with_a_deletion_inside_the_flank_records_it_as_a_flank_deletion`.

## 9. What's good

- The one-run regression (`a_one_run_bundle_scores_as_the_segment_it_is`, and its seed and
  calling counterparts) pins the bundle path to the segment's numbers bit for bit — the right
  oracle for a generalisation ([ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs)).
- The parse's lexicographic `(edits, Σ|Δcount|)` cost is additive along the path, so the
  cell-by-cell minimum is the global one — checked, and sound ([run_parse.rs](../../../../src/region_typing/run_parse.rs)).
- Slip-set enumeration is complete (the spec's own example missed `(−3, 4)`; the code does not),
  and pruning is a deterministic function of the inputs, so output stays the same on every
  platform.
- Segment selection is shared through keyed passes rather than copied, and every segment cell
  stayed identical against main on all three benchmark runs.

## 10. Commands to re-verify

All in the container (`/Users/jose/devel/pop_var_caller-str-bundles/scripts/dev.sh`):

- `cargo test --release --no-fail-fast --lib` (expect the three debug-only failures).
- `cargo test --release --lib bundle`, `… likelihood::`, `… locus_generation::`, `… psp::`.
- `bash tmp/review_fix/bench/run_all.sh <binary> <name>`, then
  `uv run --no-project python tmp/review_fix/bench/summary.py` on the host, and
  `tmp/step4/bench/segment_cells_diff.py <run dir> <name>` per run.
- HG002: `tmp/step4/hg002_two_binaries_30x.sh <out> <base> <experiment>` with
  `DEV_EXTRA_MOUNT=<main checkout>/benchmarks`, then `tmp/review_fix/region_confusion.py`.
