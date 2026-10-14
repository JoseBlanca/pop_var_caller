# Code Review: str_partial_read_nomination
**Date:** 2026-10-13
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** uncommitted change on branch `partial-candidates` (from `main` b2714614): partial reads that fit exactly one sequence of an STR locus count toward candidate selection, behind `PVC_PARTIAL_READS_NOMINATE`
**Status:** Request-changes (fixes applied the same day — `fixes_applied_str_partial_read_nomination_2026-10-13.md`)

---

## 1. Scope

- **What was reviewed:** the working-tree diff against `main` b2714614 plus the untracked
  `src/calling/allele_candidates/partial_read_votes.rs` (renamed `deciding_partial_reads.rs` by
  the fixes).
- **Reviewed against:** as-provided, worktree `pop_var_caller-partial-candidates`.
- **In-scope files:** `src/calling/allele_candidates/partial_read_votes.rs` (new),
  `src/calling/allele_candidates/ssr.rs` (`select_ssr`, `select_ssr_bundle`, `SsrSelectionConfig`),
  `src/calling/allele_candidates/mod.rs`, `src/run/callers.rs`,
  `doc/devel/ng/spec/candidate_alleles_ssr.md` §8.1, `doc/devel/ng/spec/str_bundle_loci.md` §7.6;
  and their direct callees: `leftover_of`, `compared_reads_of`, `summarise_alleles`, the record's
  `evidence_for_output`, `HeldFlank::of_witness`, the partial-read mint in
  `src/locus_generation/ssr.rs`, and the read model's censored terms in
  `src/calling/likelihood/ssr_emission.rs`.
- **Deliberately out of scope:** the segment read model's treatment of a partial read's flank bases
  (pre-existing; raised as M3 because the diff's doc claimed otherwise).
- **Categories dispatched** (owner-required deviation: no git worktrees — every reviewer read from
  this worktree with no build and no edits; and the eleven categories went to four reviewers):
  - `reliability` — always; also carried the owner's seven questions.
  - `errors`, `defaults`, `refactor_safety` — one reviewer; the switch's removal is a defaults
    question, the selection-only copy a refactor-safety one.
  - `naming`, `idiomatic`, `smells`, `module_structure` — one reviewer.
  - `float_portability`, `extras` — one reviewer; `q_sum` sums reach a written number, and the
    copy is on a per-locus path.
  - Not dispatched: `unsafe_concurrency` (no `unsafe`, no shared state), `tooling` (no
    `Cargo.toml` change).

Per-category files: `tmp/review_2026-10-13_str_partial_read_nomination/`.

## 2. Verdict

Request-changes. One output defect (B1) and the owner's adoption decision (M1) had to land before
the rule could become the default; everything else is tests, documentation and cost.

## 3. Execution status

- `cargo test --release --no-fail-fast -j 2 --lib allele_candidates` (container, before any fix)
  → exit 0, `188 passed; 0 failed`.
- `cargo fmt --check` on the four touched source files (before any fix) → exit 0.
- Clippy, the full suite and `cargo doc` were not run by the reviewers (read-only by owner
  rule); the fix run ran clippy and the touched suites.
- Findings labelled "Needs verification": 1 (Mi14).
- The reviewers ran no mutations; the fix run did (fix report §4).

## 4. Open questions and assumptions

1. **Should a counted partial read have to show at least one flank base?** Affects M2. Today a
   read that ended inside the repeat counts for the longest sequence of the table when it is
   longer than every other. That follows the rule as the owner stated it (exact fit to exactly one
   sequence of the table), and bundle reads that end inside the hull rely on it, but a sample
   carrying an expansion too long to cross would see its partial reads counted for the longest
   *spelled* sequence. Unmeasured: how many counted reads showed no flank base on the stuttering
   benchmark.
2. **Is the whole-locus edge-deletion exclusion worth refining?** Affects Mi15. A finer rule (drop
   left-held reads where some allele deletes left flank bases, right-held where some deletes right)
   is plausible but rests on two unverified assumptions; first log the held flank of the 85 reads
   behind the measured `AGC` failure, and count the loci the blunt rule silences.

## 5. Top 3 priorities

1. **B1** — `DP` counted a partial read twice whenever it counted for a sequence selection then
   refused (every refused stutter sequence at depth). Wrong numbers in every record of the adopted
   default.
2. **M1** — the switch still existed and the docs called the rule an experiment; the owner had
   decided to adopt it.
3. **B2** — no test reached a vote decided by the witness alone (every vote for a contraction) or a
   vote in a cohort of more than one sample.

## 6. Findings

### Blocker

#### B1: `src/calling/allele_candidates/ssr.rs:1337` — `FORMAT/DP` counts a deciding partial read twice when its sequence is not admitted
- **Categories:** reliability, errors, extras, smells (convergent)
- **Confidence:** High
- **Problem:** `select_ssr` and `select_ssr_bundle` rebind `observation` to the copy with the
  counted reads, and `admit_promoted_sequences_on` builds each sample's leftover
  (`UnmatchedSupport`) from it — so `num_reads` includes every partial read that counted for a
  sequence selection refused. `evidence_for_output` (`src/run/records.rs:570-583`) is handed the
  merge's own observation and adds both the leftover's `num_reads` and every partial read into
  "reads no written allele explains", which `DP` adds. On the module's own stutter fixture (32
  crossing reads, 32 deciding partial reads, the stutter sequence refused) `DP` comes out 66 for 64
  reads.
- **Fix:** take the leftover's read count and error mass from the merge's rows; ask the cap's
  "this sample earned it" test over the copy (the spec says a counted read counts wherever the bar
  is asked). Test that the leftover is the merge's.

#### B2: `partial_read_votes.rs:120` — no test reaches witness-only votes or more than one sample
- **Categories:** reliability
- **Confidence:** High
- **Problem:** every voting fixture uses an allele longer than the reference, so every voting
  partial has a whole-locus witness and a recorded side; replacing `HeldFlank::of_witness(...)` by
  `partial.held_flank` would pass the file and drop every vote for a contraction. Every fixture
  holds one sample, so zipping the votes against the wrong samples would pass too.
- **Fix:** `a_partial_read_reaching_the_flank_lifts_a_once_seen_contraction` (both sides) and
  `deciding_reads_count_only_for_the_sample_whose_reads_they_are`.

### Major

#### M1: switch still present; the docs say "not adopted"
- **Categories:** defaults, extras, naming
- **Confidence:** High
- **Problem:** `SsrSelectionConfig::partial_reads_nominate` defaults to `false`, three driver sites
  read `PVC_PARTIAL_READS_NOMINATE` (panicking on any value but 0 or 1), and §8.1, `str_bundle_loci.md`
  §7.6, the module doc and several comments call it an experiment. The owner decided on 2026-10-13
  to adopt it with no switch.
- **Fix:** delete the constant, reader, field and three call-site lines; make both call blocks
  unconditional; rewrite the tests' off arms as the same locus without its partial reads; update
  the specs and the comments that say partial reads count toward no bar (`mod.rs:1083`,
  `compared_reads_of`, `fill_sample_reads_per_rung`, `admit_promoted_sequences`' bar clause,
  `select_ssr`'s pass list). Grep showed no other `SsrSelectionConfig` built without `..at_ploidy`,
  no `PVC_` variable in `scripts/` and no caller in `tests/`; `ssr_production_differential`'s
  fixture carries no partial reads.

#### M2: `partial_read_votes.rs:170` — a read that ended inside the repeat counts for the longest sequence; the spec's prose says such a read decides nothing
- **Categories:** reliability (Major), extras (Minor)
- **Confidence:** High on the behaviour, Medium on its cost
- **Problem:** `partial_read_fits(Left)` holds whenever the read is a prefix of the allele, with no
  flank base required. With `CA×10` and `CA×12` in the table, a left-held read of `CA×11` fits only
  `CA×12` and counts for it. The precise rule (§8.1) says exactly this; the plain-words paragraph
  said the opposite. At a bundle a read ending inside the hull routinely fits one sequence only, and
  the bundle test relies on it.
- **Fix:** see open question 1. Either require a flank base (measure first) or keep it and state the
  closed reading in §8.1, pinned by a test.

#### M3: `partial_read_votes.rs:73` — at segments the vote's "fits" is not the read model's
- **Categories:** reliability
- **Confidence:** High on the code, Medium on the effect
- **Problem:** the doc said the flanks are compared "as the read model compares a partial read
  with them". At bundles they are (`flank_beyond_the_hull`, `ssr_emission.rs:2056-2065`). At
  segments `censored_from_the_tract_start` never reads a flank: it stretches the candidate to the
  read's length and compares flank bases with repeat letters, so a read that counted for `CA×12`
  because it showed `T` past it is scored as `CA×12` one unit longer with a mismatch.
- **Fix:** correct the comment and §8.1; give the segment censored term the far flank as a
  separate change.

#### M4: `partial_read_votes.rs:112` — the copy clones every partial read, which selection never reads
- **Categories:** extras (Major), idiomatic, smells (Minor)
- **Confidence:** Medium (not measured)
- **Problem:** at any locus with one counted read, the whole `CohortObservation` is cloned —
  including every sample's partial reads, one heap allocation each for their bases. Selection reads
  no partial read and no window coverage. At 1,000 samples and hundreds of reads that is tens of
  thousands of allocations per STR locus, dropped unread. It also leaves each counted read in the
  copy twice (in a row and in `partials`).
- **Fix:** copy each sample's rows with empty partials; later, reuse the rows through
  `SelectionScratch`.

#### M5: `ssr.rs:1683` — shadowing `observation` routes every later use to the copy
- **Categories:** refactor_safety
- **Confidence:** High
- **Problem:** B1 is this shadowing's first casualty; any later reader of `observation` in
  `select_ssr` gets the copy silently.
- **Fix:** name the merge's rows explicitly where they are needed.

#### M6: `partial_read_votes.rs:1` — "vote" is a third name for one concept
- **Categories:** naming
- **Confidence:** Medium
- **Problem:** the spec says a read "counts" and "singles out"; the code says "vote" and
  "deciding"; "nominate" already means a sample nominating a rung.
- **Fix:** name the module and functions after the spec's words.

### Minor

#### Mi1: `src/run/explain.rs:231` — `--explain-loci` shows rows that do not explain the selection
- **Categories:** reliability, errors · **Confidence:** High
- The explain rows are built from the merge's rows for a selection made over the copy, so an allele
  admitted on one crossing read and one counted partial read shows one read under a floor of two.
  Fix: an explain row per counted read.

#### Mi2: `partial_read_votes.rs:92` — locus length derived differently from the read model's
- **Categories:** reliability, errors · **Confidence:** High
- `alleles[0].len()` against `LocusLen::of_region(region)`; equal today, and a divergence would read
  the read from the wrong side silently. Fix: `LocusLen::of_region`.

#### Mi3: `partial_read_votes.rs:197` — the flank-deletion slicing is unreachable and untested
- **Categories:** reliability, refactor_safety · **Confidence:** High
- The edge-deletion guard returns first, so `flank_deletion` is always `NONE` there. Fix: test
  `partial_read_fits` directly with a left and a right deletion.

#### Mi4: `partial_read_votes.rs:133` — the error mass a counted read adds is never asserted
- **Categories:** reliability · **Confidence:** High

#### Mi5: `partial_read_votes.rs:93` — a sample with no crossing read can now reach the bar on two partial reads
- **Categories:** reliability, errors (cross-category) · **Confidence:** High
- Allowed by the precise rule, unstated in the spec, untested. Fix: state and pin it.

#### Mi6: `partial_read_votes.rs:85` — nothing records how many reads counted, or how many loci the edge-deletion guard silenced
- **Categories:** errors · **Confidence:** High

#### Mi7: `partial_read_votes.rs:141` — `..AlleleSupport::default()` in a row a counted read creates
- **Categories:** errors, refactor_safety · **Confidence:** High
- Strand count, both mapping-quality sums and `placed_left` are zero; nothing reads them from the
  copy today, but a new `AlleleSupport` field would be zero-filled with no compile error.

#### Mi8: `ssr.rs:2406` and others — "partial reads count nowhere" tests pass only because `partial_of`'s bases fit no allele
- **Categories:** refactor_safety · **Confidence:** High

#### Mi9: stale docs
- **Categories:** reliability, naming, defaults · **Confidence:** High
- `PartialObservation::held_flank` (`build.rs:1498-1502`) and `Classified::Observed::held_flank`
  (`locus_generation/ssr.rs:942-944`) say the side is recorded only at bundles; the mint records it
  at segments too (`ssr.rs:1753-1754`).

#### Mi10: function and local names (`with_deciding_partials_counted`, `the_one_allele_a_partial_read_fits`, `far`, `past`, `first`, `held`)
- **Categories:** naming · **Confidence:** Medium

#### Mi11: `partial_read_votes.rs:197` — the far-flank slicing duplicates `flank_beyond_the_hull`
- **Categories:** smells · **Confidence:** High

#### Mi12: `partial_read_votes.rs:75` — one 77-line function does three jobs
- **Categories:** smells · **Confidence:** Medium

#### Mi13: the discovery pre-pass also reads the counted rows; §8.1 did not say so
- **Categories:** errors (cross-category) · **Confidence:** High

#### Mi14: mate overlap — Needs verification
- **Categories:** reliability · **Confidence:** Low
- No mate-overlap deduplication was found in the STR generator, so a fragment whose one mate
  crosses and whose other is a deciding partial read now counts twice where it counted once. Two
  crossing mates already counted twice, so this extends a convention rather than creating one.

#### Mi15: the edge-deletion exclusion is whole-locus
- **Categories:** reliability · **Confidence:** Medium. See open question 2.

#### Mi16: `mod.rs:34` — the module is a public peer of the two selection paths
- **Categories:** module_structure, idiomatic · **Confidence:** Medium
- It is one step of the STR path; nothing outside `allele_candidates` names it.

### Nits

- The `read.len() <= allele.len() + far.len()` guard is redundant (`starts_with`/`ends_with` with a
  longer needle already fail).
- The `sequence_may_run_on_past` guard is dead at STR loci (the re-read is the generic path's);
  say so.
- `partial_reads_nominate_from_environment` has no test of its panicking arms (moot after M1).
- The edge-deletion test covers only a left deletion.
- `Vec<Vec<…>>` per locus where a flat list would do; `witnessed_in_locus.clone()` per partial
  (a `SmallVec` copy, not an allocation).
- Tests use `super::super::fixtures`.
- `examples/ng_candidate_selection_probe.rs` changes output once the rule is unconditional; figures
  measured with it before are stale.

## 7. Out of scope observations

- The segment read model ignores a partial read's flank bases (`censored_from_the_tract_start`) —
  a likelihood change of its own (M3).

## 8. Missing tests to add now

Grouped by function; all written in the fix run.

- `count_deciding_partial_reads`: `a_partial_read_reaching_the_flank_lifts_a_once_seen_contraction`
  (B2), `deciding_reads_count_only_for_the_sample_whose_reads_they_are` (B2),
  `a_read_ending_inside_the_repeat_past_every_shorter_allele_counts_for_the_longest` (M2),
  `a_sample_without_crossing_reads_nominates_on_two_deciding_partial_reads` (Mi5), the right-flank
  case in `no_partial_read_counts_where_an_allele_deletes_across_the_edge` (Nit), the error-mass
  asserts in `deciding_partials_leave_a_stutter_sequences_share_where_it_was` (Mi4).
- `select_ssr` / `leftover_of`: `the_leftover_holds_the_merges_reads_and_not_the_counted_partials` (B1).
- `partial_read_fits`: `partial_read_fits_takes_the_alleles_own_flank_deletion_from_the_far_side`
  (Mi3), `partial_read_fits_at_its_boundaries`.
- Not written: a property test that the right-held arm mirrors the left-held one.

## 9. What's good

- The deciding reads enter numerator and denominator through the same rows, so every place the bar
  is asked (`summarise_alleles`, rung nomination, discovery, the cap's earned test) sees one
  consistent count — `deciding_partials_leave_a_stutter_sequences_share_where_it_was` pins it.
- The `CohortObservation` copy is an exhaustive struct literal: a new field breaks the build there
  rather than being silently dropped.
- The edge-deletion guard is documented with the measured failure that motivated it (one `AGC`
  segment, 18 right calls lost).
- `counted_rows_keep_the_merges_order` pins the ordering invariant every walk over the rows asserts.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --no-fail-fast -j 2 --lib allele_candidates`
- `scripts/dev.sh cargo test --release --no-fail-fast -j 2 --lib cross_platform_digests`
- The fix run's identity check: `DEV_EXTRA_MOUNT=<str-bundles>/tmp/bundle_bench/runs scripts/dev.sh tmp/verify/run.sh`
