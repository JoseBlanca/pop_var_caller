# Code Review: str_bundles_simulator
**Date:** 2026-10-11
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** uncommitted change on branch `str-bundles` (against ccef7dca) that makes the stuttering-benchmark simulator lay STR bundles, and makes the driver, the typing check and the scorer handle them
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the working-tree diff against HEAD ccef7dca, plus the one new untracked file.
- **Reviewed against:** branch `str-bundles`, HEAD ccef7dca, uncommitted changes as of 2026-10-11.
- **In-scope files:**
  - `examples/ng_str_stutter_benchmark.rs` (the simulator; +1,276 lines)
  - `examples/ng_typed_region_dump.rs` (the typing dump; gains a `members` column)
  - `scripts/ng_str_stutter_check_bundles.py` (new: the bundle typing check)
  - `scripts/ng_str_stutter_score.py` (bundle cells, a third vcfeval stretch)
  - `scripts/ng_str_stutter_benchmark.sh` (runs the bundle check; HipSTR left off bundles)
  - `scripts/ng_str_accuracy/copy_run.sh` (copies the bundle files)
  - `doc/devel/ng/spec/str_stutter_benchmark.md` §11 (new), `doc/devel/ng/spec/str_bundle_loci.md` §9
- **Out of scope:** `src/` (untouched by the change; read only to check the scanner copy against
  `src/tandem_repeat.rs` and `src/region_typing/segment_criteria.rs`). Readers downstream of the
  scorer (`scripts/ng_str_accuracy/compare.py`, `score_callers.py`, `edge_deletions.py`,
  `scripts/ng_str_stutter_grid_summary.py`, `scripts/ng_str_stutter_true_parameters.py`) were read
  to confirm they are unaffected, not reviewed.
- **Deviation from the skill, required by the owner's disk budget:** no git worktrees were made.
  Two read-only sub-agents reviewed from this worktree without building or editing; the
  orchestrator did every build and every run.
- **Categories dispatched:**
  - `reliability`, `extras` (diff matches stated intent) — orchestrator, with builds and runs:
    this is a fixture whose correctness is in its coordinates.
  - `naming`, `idiomatic`, `smells`, `errors`, `refactor_safety`, `module_structure`,
    `float_portability`, `defaults` — one read-only sub-agent over the two Rust files
    (`tmp/review_2026-10-11_str-bundles-sim/rust_categories.md`).
  - The Python and shell scripts, all categories — one read-only sub-agent
    (`tmp/review_2026-10-11_str-bundles-sim/scripts_review.md`).
  - `unsafe_concurrency`, `tooling` — skipped: the change adds no shared state, threads or
    dependencies.

## 2. Verdict

**Approve-with-changes.** Every one of the implementer's six claims holds when measured (§3). No
finding changes a number the benchmark has produced. What needs fixing is a stale-cache hazard in
the scorer that can count bundle-window variants twice, and test cover for the coordinate arithmetic
the bundle alleles depend on, which today is correct but unprotected.

## 3. Execution status

Commands run (all inside the dev container):

| command | exit | result |
|---|---|---|
| `cargo build --release --example ng_str_stutter_benchmark --example ng_typed_region_dump --bin pvc`, plus HEAD's simulator built as a temporary example (`zz_review_before_bundles`, deleted afterwards) | 0 | built |
| `cargo test --release --example ng_str_stutter_benchmark` | 0 | `test result: ok. 14 passed; 0 failed` |
| `tmp/review_2026-10-11_str-bundles-sim/verify.sh` (output saved beside it as `verify_before_fixes.txt`) | 0 | below |

What `verify.sh` measured, claim by claim:

1. **`bundles=0` is byte-identical to the build before the change — CHECKED-CORRECT** on a
   configuration unlike the implementer's (`design=heavy samples=3 depth=2 segments_per_period=60
   genome_bp=500000 contigs=3 inbreeding=0.1 seed=7`): all 14 files the old build writes have the
   same SHA-256, the six FASTQs compared decompressed. The new build adds the three bundle files,
   header only. The simulator's standard output differs (a `STR bundles` line, and `bundles: 0` in
   the settings line), and the driver saves it as `simulate.log` (see Mi5).
2. **The scanner copy rejects shapes the typing would not see as drawn — CHECKED-CORRECT** on two
   seeds other than the default: seed 1 and seed 424242 each typed 500 of 500 bundles as one bundle
   with exactly the drawn members (tolerance 0, largest offset 0), no bundle the simulator did not
   lay, and all 600 segments as planted.
3. **Gap and hull rules — CHECKED-CORRECT** by reading `draw_gap` and the shape test, and by the
   typing check above. Flanks are 150 bases either side, spacer included.
4. **One stutter site per member at its own rate — CHECKED-CORRECT** by reading `build_haplotype`,
   `cut_molecule` and `draw_bundle_alleles`: each allele computes its member runs from its own
   counts, so a longer first member moves the later runs; the left edge deletion shifts every later
   run by all the deleted bases, the right one shortens only the last run. Not covered by a test
   (M2).
5. **Separate files; segment outputs unchanged; `outside` excludes bundle windows —
   CHECKED-CORRECT.** The scripts reviewer scored seed 1 with HEAD's scorer and the new one:
   `score_cells.tsv` and `score_summary.tsv` byte-identical. The three vcfeval stretches cover every
   reference base exactly once on seed 1. Scoring the truth VCF as a caller (`SCORE_SELF_CHECK=1`)
   gives `right` or `ref_right` in every bundle cell: seed 1, 632 right and 868 ref_right of 1,500;
   seed 424242, 617 and 883.
6. **HipSTR left off bundles — CHECKED-CORRECT** (region file from `segments.bed` only; the scorer
   skips HipSTR for bundle cells and the bundle stretch).

The distribution figures in `str_stutter_benchmark.md` §11.1, recomputed at the default seed: 81
three-member bundles of 500; members by motif length 248, 185, 235, 167, 143, 103; 93 of 581
neighbouring pairs share a motif length; hull median 77, 90th percentile 96 — all CHECKED-CORRECT.
The explanation given for the median is incomplete (Mi4).

Commands not run: `cargo fmt --check`, `cargo clippy --all-targets`, `cargo doc`, `cargo audit`, the
full test suite — the owner's light dev loop (touched tests only); the change touches no library
code. vcfeval was not run (no RTG in the image at the run path; the scorer was run with
`SKIP_VCFEVAL=1`).

Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. **Should the bundles' shape distribution follow HG002's?** (affects Mi4.) The simulator's median
   hull is 77 bases against HG002's 35, because each member's count is drawn uniform from the
   calling floor to the segments' highest (30 copies for motif lengths one and two, 15 for the
   rest), so the average member is 31 bases. A long hull leaves fewer reads that cross it, so at 3–5×
   the benchmark will be harder on bundles than HG002 is. The review only corrects the explanation;
   whether to change the draw is the owner's call.
2. **Is reusing an older run's vcfeval results a case to support?** (affects M1.) The fix assumes yes
   and invalidates by input contents.

## 5. Top 3 priorities

1. **M1** — the scorer's vcfeval cache can reuse an `outside` result computed before bundle windows
   were excluded, counting each bundle-window variant twice with no warning.
2. **M2** — no test pins where a bundle allele's member runs and stutter stretches sit; the edge
   deletion and the shift of later members are correct today and unguarded.
3. **M3** — the scanner copy takes its two scores from the library while its doc says the fixture
   does not follow the library; and nothing but a stopped benchmark run reports the copy drifting.

## 6. Findings

### Major

#### M1: `scripts/ng_str_stutter_score.py:412-417` — vcfeval results are cached by stretch name, so a stale `outside` result is reused
- **Confidence:** High for the mechanism; Medium that a run in that state exists.
- **Categories:** scripts review (orchestrator confirmed by reading).
- **Problem:** `run_vcfeval` skips a stretch whenever `vcfeval/<caller>.<sample>.<stretch>/summary.txt`
  exists, whatever the BED, truth or call file now holds. Input: a run simulated with bundles and
  scored once by the scorer at HEAD (its `outside` covers the bundle windows), then rescored with
  `REDO=score`. Output: the old `outside` counts are reused, a new `bundles` stretch is computed,
  and `score_vcfeval.tsv` counts every variant in a bundle window twice. The same lookup keeps an old
  `truth.vcf.gz` after `REDO=simulate` (rebuilt only when its index is missing).
- **Why it matters:** wrong SNP/indel counts outside the STR windows, with nothing in the output to
  say so.
- **Suggested fix:** store a digest of each stretch's inputs (the BED, the truth, the call file)
  beside its `summary.txt` and redo the stretch when the digest differs; rebuild `truth.vcf.gz`
  when `truth.vcf` is newer.

#### M2: `examples/ng_str_stutter_benchmark.rs:1395-1568` — the bundle alleles' member runs and stutter stretches are untested
- **Confidence:** High.
- **Categories:** reliability.
- **Problem:** `bundle_alleles_spell_their_member_counts` checks member runs for pure length alleles
  only; `a_bundle_has_one_stutter_site_per_member` checks each site's motif length and count but not
  where it lies. The coordinates that move — later members after a length change, every member
  after a left edge deletion (`run.0 - deleted`), the last member's end after a right one — are
  unchecked. A wrong offset there (for example `run.0 - in_repeat` at line 1541) puts stutter in a
  gap or a flank: the reads change, no test fails, and the truth VCF stays right, so even the
  scorer's self-check passes.
- **Why it matters:** the benchmark exists to measure stutter; misplaced stutter is a wrong
  fixture that looks right.
- **Suggested fix:** a test over every allele kind: between consecutive runs lie exactly the gap
  bases; an untouched run is its motif `count` times; the cut member is a periodic stretch holding
  `count` whole copies and ending (left cut) or starting (right cut) on a whole motif; a
  substitution changes one base; each pure span lies inside its run and repeats with the member's
  motif length.

#### M3: `examples/ng_str_stutter_benchmark.rs:123, 833-891` — the scanner copy follows two library constants and has no test against the library
- **Confidence:** High.
- **Categories:** refactor_safety (Rust reviewer), reliability (orchestrator).
- **Problem:** `scanner_repeats`' doc says it is the file's own copy "so the fixture does not change
  when the caller's code does", but its scores are `pop_var_caller::tandem_repeat::{DEFAULT_MATCH_REWARD,
  DEFAULT_MISMATCH_PENALTY}`. Retuning either in `src/tandem_repeat.rs` changes which shapes are
  accepted, so every bundle and every draw after them, with no test failing (the hash test runs at
  `bundles=0`). In the other direction, if the library's arithmetic changes (for example its
  tie-breaking between equal-scoring stretches, which the copy reimplements recursively rather than
  with Ruzzo–Tompa's stack), the copy keeps accepting shapes the typing now splits differently,
  and the first sign is the driver stopping at its typing check.
- **Why it matters:** claim 2's guarantee rests on the copy and the library agreeing, and nothing
  checks that they do except a full benchmark run.
- **Suggested fix:** pin the two scores and the gap bound locally, with a test asserting they equal
  the library's; and a differential test running the copy and `find_tandem_repeats` (plus the copy
  floor) on a few thousand bundle-like sequences whose gaps follow no rule.

### Minor

#### Mi1: `scripts/ng_str_stutter_score.py:378` — `outside` overlaps a window when two windows nest
- **Confidence:** High (reproduced by the scripts reviewer: segment window 100–400 with bundle window
  200–300 gives `outside` 300–1000, so bases 300–400 count as outside).
- **Problem:** `cursor = end` rather than `max(cursor, end)`. Only latent: segment and bundle windows
  are at least 843 bases apart on both seeds.
- **Suggested fix:** `cursor = max(cursor, end)`.

#### Mi2: `examples/ng_str_stutter_benchmark.rs:3060` — the truth-VCF test's segment counter counts bundle records
- **Confidence:** High.
- **Problem:** `segment_records += usize::from(fields[2] != ".")` counts `bun…` IDs too, so
  `assert!(segment_records > 0)` passes with no segment record.
- **Suggested fix:** count `starts_with("seg")` and `starts_with("bun")` separately.

#### Mi3: `examples/ng_str_stutter_benchmark.rs:41-43` — "byte for byte" reads as holding at every setting
- **Confidence:** High.
- **Problem:** the segments' files keep their columns and hold segments only; they are what they
  were byte for byte only at `bundles=0`. At the default 500 the segments sit elsewhere.
- **Suggested fix:** say both halves.

#### Mi4: `doc/devel/ng/spec/str_stutter_benchmark.md` §11.1 — the hull median's cause is misstated
- **Confidence:** High (measured: 1,081 members at the default seed average 30.9 bases).
- **Problem:** "HG002's is 35, because members here reach 30 copies". Only motif lengths one and two
  reach 30 copies; the cause is that every member's count is uniform up to the segments' highest,
  so long-motif members of 15 copies (45–90 bases) count too.
- **Suggested fix:** state the mechanism and the average member length; leave open question 1 to
  the owner.

#### Mi5: `scripts/ng_str_stutter_first_pass.sh`, `simulate.log` — the new default changes what old commands produce
- **Confidence:** High.
- **Problem:** `bundles` defaults to 500, so rerunning the first-pass grid (which passes no
  `bundles=`) produces cohorts other than the 2026-10-10 report's: the segments move. And
  `FIT_FROM` compares `simulate.log` settings lines, which now carry `bundles:`, so a fit from a
  run made before bundles is refused even at `bundles=0`. The refusal is safe; the silent cohort
  change is not visible anywhere a reader of the grid script would look.
- **Suggested fix:** say in the grid script's header that runs since 2026-10-11 lay 500 bundles and
  that the 2026-10-10 cohorts need `bundles=0`.

#### Mi6: `scripts/ng_str_stutter_score.py:433-456` — `longest_period` reads as "the largest motif length"
- **Confidence:** High (seed 1, `bun00000` = ATCT×13 + AAAAC×6: `longest_period` 4, largest motif
  length 5).
- **Suggested fix:** name the columns `longest_member_period` and `longest_member_repeat_band`.

#### Mi7: `examples/ng_str_stutter_benchmark.rs:1418-1432, 1584` — segment and bundle allele draws duplicate their constants
- **Confidence:** High. **Categories:** smells.
- **Problem:** the step loop, `2 + below(5)`, `0..50`, `2 + below(9)`, the substitution clamp and
  `0.2 + 0.6 * unit` are written twice; a change to one copy gives segments and bundles different
  populations silently.
- **Suggested fix:** shared helpers and named constants, keeping the draw order (the hash test pins it).

#### Mi8: `examples/ng_str_stutter_benchmark.rs:1369-1382` — per-member facts in parallel vectors joined by `zip`
- **Confidence:** Medium. **Categories:** smells, idiomatic.
- **Problem:** `repeats`, `member_runs`, `pure_spans` and `shape.periods()` are zipped; a short
  vector silently drops stutter sites. `pure_spans` is empty until a second pass fills it.
- **Suggested fix:** one `Vec<AlleleMember { repeats, run, pure_span }>`.

#### Mi9: `examples/ng_str_stutter_benchmark.rs:593, 846` — members and scanner hits are anonymous tuples
- **Confidence:** High. **Categories:** naming, idiomatic.
- **Suggested fix:** `BundleMember { motif, reference_repeats }` and `ScannedRepeat { start, end, period }`.

#### Mi10: `examples/ng_typed_region_dump.rs:151-166` — a second copy of `members_json`, behind a wildcard arm
- **Confidence:** High. **Categories:** module_structure, refactor_safety.
- **Problem:** the format duplicates the private `src/cli/typed_regions.rs::members_json`; and
  `_ => "."` over `RegionKind` would absorb a new variant silently.
- **Suggested fix:** spell the three other arms; making one public formatter shared by both is a
  library API change for the owner.

#### Mi11: `examples/ng_str_stutter_benchmark.rs:734` — the shape redraw has no bound
- **Confidence:** Medium. **Categories:** errors.
- **Problem:** if a later constant change makes a three-member shape unsatisfiable, the simulator
  spins forever with no message, where its sibling helpers panic with one.
- **Suggested fix:** a bounded loop that panics naming the rule; the accepted shapes are unchanged.

#### Mi12: `scripts/ng_str_stutter_score.py:436-440` — each caller's VCF is parsed twice
- **Confidence:** High.
- **Problem:** `score_bundles` re-reads every call file the segment loop already read. Runtime only.
- **Suggested fix:** score bundles inside the segment loop's caller iteration.

#### Mi13: `scripts/ng_str_stutter_score.py` docstring — `outside` changed meaning under its old name
- **Confidence:** High.
- **Problem:** in a run with bundles `outside` excludes the bundle windows; in an older run it does
  not, and `grid_vcfeval.tsv` mixes the two under one name.
- **Suggested fix:** keep the name, which readers filter on, and say in the docstring that it
  differs between runs with and without bundles.

### Nits

- `scripts/ng_str_stutter_check_bundles.py:53`: a dump built before the `members` column fails with a
  bare `IndexError`; members are compared in dump order rather than sorted by start; the overlap
  scan reaches back only 1,000 bases.
- `scripts/ng_str_stutter_score.py:499`: a run whose bundle table has a header and no rows leaves any
  older `score_bundle_*.tsv` in place.
- `scripts/ng_str_accuracy/score_callers.py:7-8`: its docstring does not mention the bundle files.
- `examples/ng_str_stutter_benchmark.rs`: `shapes`/`segment_bp` now cover bundles; `joined` is
  shadowed by a test local; `[2; 6]` (the clean-flank cap) appears four times; `first`/`second`/
  `moved` closures; `fnv` → `fnv1a_64`; "one for every six segments" holds at the default count only;
  the realised share of two-member changes is a little under a quarter, because a second move below
  one copy is dropped.
- The probe in `scanner_reports_exactly_the_members` lays its left spacer on an empty prefix, so it
  can differ from the spacer laid on the contig; harmless, since any spacer base breaks the member,
  and the typing check is the backstop.

## 7. Out of scope observations

- `scripts/ng_str_stutter_grid_summary.py` reads no bundle files, so bundle results reach no grid
  table yet. A follow-up when bundle runs are summarised.
- `PROJECT_STATUS.md` has no block for the STR-bundles work; its current-focus section stops at
  2026-10-10.

## 8. Missing tests to add now

- `every_bundle_allele_keeps_its_member_runs_on_its_members` (M2) — every allele kind; catches a
  wrong shift of later members or a wrong edge-deletion run.
- `the_scanner_copy_agrees_with_the_library_scanner` (M3) — 3,000 bundle-like sequences with
  unruled gaps, copy against `find_tandem_repeats` plus the copy floor; catches drift in either.
- `the_scanner_copy_uses_the_library_scores` (M3) — the pinned scores and gap bound against the
  library's constants.
- In `truth_vcf_reproduces_the_truth_sequences`, separate segment and bundle record counts (Mi2).

## 9. What's good

- The bundles' randomness lives in a second stream (`draw_cohort`), which is what makes
  `bundles=0` byte-identical, and a hash test pins it rather than a prose promise.
- `WindowRecord` (`examples/ng_str_stutter_benchmark.rs:2123`) writes segment and bundle records
  through one trimming and phasing path, so the truth-VCF test covers both.
- The typing check (`scripts/ng_str_stutter_check_bundles.py`) compares members at tolerance 0 and
  also reports typed bundles nobody laid — it fails in both directions.
- Edge-deletion coordinates carry a comment saying which runs move and by how much (line 1535).

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --example ng_str_stutter_benchmark`
- `scripts/dev.sh bash tmp/review_2026-10-11_str-bundles-sim/verify.sh` (expects the binaries in
  `tmp/review_2026-10-11_str-bundles-sim/bins/`, built by `build.sh` beside it)
