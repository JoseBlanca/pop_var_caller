# Fix Application Report: str_partial_read_nomination_2026-10-13.md

**Date:** 2026-10-13
**Source review:** `doc/devel/reports/reviews/str_partial_read_nomination_2026-10-13.md`
**Source state reviewed against:** uncommitted change on branch `partial-candidates`, from `main` b2714614
**Execution mode:** non-interactive (owner's decision of 2026-10-13 given up front: adopt the rule as the default and remove the switch)
**Overall status:** Completed (uncommitted, as instructed)

---

## 1. Executive summary

The rule is now the caller's only behaviour: a partial read whose bases fit exactly one sequence
of an STR segment's or bundle's table counts as a read of it, in both halves of the "two reads and
a tenth" bar. There is no environment variable and no flag. The one output defect the review found
is fixed: the record's `DP` counted such a read twice when the sequence it counted for was refused.
Re-called with the adopted default, the two benchmark runs reproduce the earlier switch-on VCFs in
every genotype, allele and `AD`; only `DP` moves, and only down (§4, B1). The measured accuracy
numbers in `candidate_alleles_ssr.md` §8.1 therefore stand, and no run was re-measured.

### Review totals
- Blockers: 2
- Majors: 6
- Minors: 16
- Nits: 7 (grouped)

### Outcome totals
- Applied: 13 (B1, B2, M1, M4, M6, Mi2, Mi3, Mi4, Mi5, Mi9, Mi10, Mi12, Mi13)
- Applied with adaptation: 3 (M2, M3, M5)
- Already fixed: 0
- Deferred: 7 (Mi1, Mi6, Mi7, Mi8, Mi11, Mi14, Mi15)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 1 (Mi16, by M6)
- Awaiting user answer: 0 (two open questions recorded for the owner, below)

### Validation summary
All in the dev container, release profile, `-j 2`:
- `cargo fmt --check` on the seven touched source files → exit 0.
- `cargo clippy --release -j 2 --lib --tests` → no warning or error in any touched file (the
  vendored `noodles-cram` warnings are pre-existing).
- `cargo test --release --no-fail-fast -j 2 --lib allele_candidates` → `195 passed; 0 failed`
  (188 before the fixes, 7 of them in this module; 7 tests added).
- `… --lib run::` → `779 passed; 0 failed`.
- `… --lib calling::` → `1248 passed; 0 failed; 3 ignored`.
- `… --lib cross_platform_digests` → `4 passed; 0 failed`. **The checksums did not move**, so
  nothing was re-recorded.
- The full `--all-targets` suite, `cargo doc` and `cargo audit` were not run: touched-tests-only is
  the owner's dev-loop rule, and a full test build OOM-killed a 16 GB container before.
- Performance check: skipped — nothing in `benches/` reaches candidate selection.

### Unresolved high-priority findings
- None of Blocker or Major severity. Two owner questions remain (§3).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | `DP` counts a deciding partial read twice | Apply | Applied | No | `mod.rs`, `ssr.rs`, `generic.rs`, `deciding_partial_reads.rs` | Pass; mutant killed; benchmark DP check | No |
| B2 | Blocker | No test reaches witness-only votes or several samples | Apply | Applied | No | `deciding_partial_reads.rs` | Pass; 2 mutants killed | No |
| M1 | Major | Switch still present; docs say not adopted | Apply | Applied | Owner decision | `ssr.rs`, `callers.rs`, `mod.rs`, module, both specs | Pass | No |
| M2 | Major | Read ending inside the repeat counts for the longest sequence | Apply (docs + test) | Applied with adaptation | No | spec §8.1, module doc, test | Pass | Yes (question 1) |
| M3 | Major | Segment read model compares partial reads differently | Apply (docs) | Applied with adaptation | No | module doc, spec §8.1 | N/A | Yes (likelihood change) |
| M4 | Major | Copy clones every partial read | Apply | Applied | No | `deciding_partial_reads.rs` | Pass; asserted | Scratch reuse later |
| M5 | Major | Shadowed `observation` reaches the leftover | Apply | Applied with adaptation | No | `ssr.rs` | Pass | No |
| M6 | Major | "Vote" is a third name | Apply | Applied | No | module renamed, functions renamed | Pass | No |
| Mi1 | Minor | `--explain-loci` rows do not explain the selection | Defer | Deferred | No | None | N/A | Yes |
| Mi2 | Minor | Locus length derived differently | Apply | Applied | No | `deciding_partial_reads.rs` | Pass | No |
| Mi3 | Minor | Flank-deletion slicing untested | Apply | Applied | No | test | Pass; mutant killed | No |
| Mi4 | Minor | Error mass never asserted | Apply | Applied | No | test | Pass | No |
| Mi5 | Minor | A sample with no crossing read can reach the bar | Apply | Applied | No | spec, test | Pass | No |
| Mi6 | Minor | No count of what the rule did | Defer | Deferred | No | None | N/A | Yes |
| Mi7 | Minor | `..AlleleSupport::default()` in a created row | Defer | Deferred | No | comment only | N/A | Yes |
| Mi8 | Minor | Older "partials count nowhere" tests pass by fixture accident | Defer | Deferred | No | None | N/A | Yes |
| Mi9 | Minor | Stale docs | Apply | Applied | No | `build.rs`, `locus_generation/ssr.rs`, `mod.rs`, `ssr.rs` | fmt | No |
| Mi10 | Minor | Function and local names | Apply | Applied | No | module | Pass | No |
| Mi11 | Minor | Far-flank slicing duplicates `flank_beyond_the_hull` | Defer | Deferred | No | None | N/A | Yes |
| Mi12 | Minor | One function does three jobs | Apply | Applied | No | module | Pass | No |
| Mi13 | Minor | Discovery pre-pass reads the counted rows, unstated | Apply | Applied | No | spec §8.1 | N/A | No |
| Mi14 | Minor | Mate overlap may count a fragment twice | Defer | Deferred | No | None | N/A | Yes |
| Mi15 | Minor | Edge-deletion exclusion is whole-locus | Defer | Deferred | No | None | N/A | Yes (question 2) |
| Mi16 | Minor | Module is a public peer of the selection paths | Apply | Superseded (by M6) | No | `mod.rs` | Pass | No |
| Nits | Nit | grouped | Apply most | Applied (5 of 7) | No | module | Pass | No |

## 3. Questions asked and answers

None asked during the run; the owner's adoption decision was given before it. Two questions are
left for the owner, carried from the review:

1. **M2** — Should a counted partial read have to show at least one base of the flank past the
   allele? Today it need not, which follows the rule as stated and lets bundle reads that end inside
   the hull count; the risk is a sample with an expansion too long to cross. Unmeasured.
2. **Mi15** — Refine the edge-deletion exclusion to the side the deletion is on? Plausible, not
   clearly safe; needs the held flank of the 85 reads behind the measured `AGC` failure first.

## 4. Per-finding log

### B1 — `DP` counts a deciding partial read twice
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** the selection copy's rows reached `UnmatchedSupport::num_reads`, which the record
  adds to `DP` beside every partial read.
- **Implementation summary:** `leftover_of` takes two row sets: the merge's (`sample`), from which
  the read count and error mass come, and the rows selection read (`selection_rows`), over which
  the cap's "this sample earned it" test asks the bar. `admit_promoted_sequences{,_on}` take
  `merge_rows: &[SampleSupport]`, which `select_ssr` and `select_ssr_bundle` pass from the merge's
  observation; the generic path passes the same rows twice. Field and function docs say which
  count includes counted reads.
- **Review suggestion used verbatim?:** No (same design, written against the code).
- **Adaptation:** None of substance.
- **Verification performed:**
  - test-first by mutation: with the leftover taken from the copy again, exactly
    `the_leftover_holds_the_merges_reads_and_not_the_counted_partials` fails (194 pass, 1 fail);
  - on the 50-sample, PCR, 5× benchmark run: against the earlier switch-on VCF, every genotype,
    allele and `AD` is identical; `INFO/DP` is lower in 1,154 of 7,829 records and `FORMAT/DP` in
    1,653 sample cells, by 1,668 reads in all, and higher nowhere. At one sample, PCR, 5× the VCF
    body is byte-identical.
- **Files changed:** `src/calling/allele_candidates/{mod.rs,ssr.rs,generic.rs,deciding_partial_reads.rs}`
- **Tests added or modified:** `the_leftover_holds_the_merges_reads_and_not_the_counted_partials`
- **Validation:** `cargo test --release -j 2 --lib allele_candidates` → 0, 195 passed.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** `earned_reads_cut_by_the_cap` may now exceed `num_reads` (it counts counted
  partial reads); only whether it is zero is read, and the field doc says so.

### B2 — No test reaches witness-only votes or several samples
- **Severity:** Blocker · **Final status:** Applied
- **Implementation summary:** `a_partial_read_reaching_the_flank_lifts_a_once_seen_contraction`
  (left- and right-held reads shorter than the locus, asserting their side is in the witness alone)
  and `deciding_reads_count_only_for_the_sample_whose_reads_they_are` (two samples).
- **Verification performed:** mutants — side read from the recorded flank only: only the
  contraction test fails; votes zipped against reversed samples: only the two-sample test fails.
- **Validation:** as B1.

### M1 — Switch still present
- **Severity:** Major · **Final status:** Applied · **User input:** the owner's decision of 2026-10-13
- **Implementation summary:** deleted `PARTIAL_READS_NOMINATE`, `partial_reads_nominate_from_environment`,
  `SsrSelectionConfig::partial_reads_nominate` and its default, and the three driver lines in
  `callers.rs` (that file is now back to `main`); both call blocks unconditional. Tests' off arms
  became "the same locus without its partial reads". Docs that said partial reads count toward no
  bar now say which do: `compared_reads_of`, `summarise_alleles`' empty-sample comment,
  `fill_sample_reads_per_rung`, `admit_promoted_sequences`' bar clause, `select_ssr`'s pass list.
  `candidate_alleles_ssr.md` §8 and §8.1 rewritten as adopted (dated note, measured table kept,
  reproduction recorded); `str_bundle_loci.md` §7.2, §7.6 and §11 item 5 updated.
- **Verification performed:** `grep` for the variable and the field across `src/`, `examples/`,
  `tests/`, `scripts/`, `doc/devel/ng/` finds nothing.
- **Follow-up:** `examples/ng_candidate_selection_probe.rs` now prints with the rule on; figures it
  printed before are of the old selection.

### M2 — Read ending inside the repeat counts for the longest sequence
- **Severity:** Major · **Final status:** Applied with adaptation
- **Reasoning:** the behaviour is the rule as the owner stated it (exact fit to exactly one sequence
  of the table); bundle reads that end inside the hull rely on it, so requiring a flank base would
  change measured behaviour on a question the owner has not been asked. Changing it would also have
  required re-measuring every run.
- **Implementation summary:** behaviour kept; the spec's and module's prose now state the closed
  reading and its risk; pinned by
  `a_read_ending_inside_the_repeat_past_every_shorter_allele_counts_for_the_longest`.
- **Follow-up:** open question 1.

### M3 — Segment read model compares partial reads differently
- **Severity:** Major · **Final status:** Applied with adaptation (documentation)
- **Implementation summary:** the false "as the read model compares" claim is removed; §8.1 gains a
  paragraph saying the bundle model reads the far flank and the segment model does not, and what
  that does to a counted read's weight. The likelihood change is out of this diff's scope.
- **Follow-up:** give `censored_from_the_tract_start` the far flank, as a separate change.

### M4 — Copy clones every partial read
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** the copy is built sample by sample with an exhaustive destructure,
  cloning the support rows and carrying `partials: Vec::new()`; asserted in
  `deciding_partials_leave_a_stutter_sequences_share_where_it_was`.
- **Follow-up:** reusing the rows through `SelectionScratch` is not done; cost unmeasured.

### M5 — Shadowed `observation`
- **Severity:** Major · **Final status:** Applied with adaptation
- **Adaptation:** the shadow is kept — every pass after it is meant to read the copy — and the one
  consumer that must not, the leftover, now receives `merge_rows` by name, bound before the shadow.

### M6 — "Vote" is a third name
- **Severity:** Major · **Final status:** Applied
- **Implementation summary:** module renamed `partial_read_votes` → `deciding_partial_reads` and
  made private (`mod`, not `pub mod`); `with_deciding_partials_counted` → `count_deciding_partial_reads`;
  the row-adding step extracted as `add_deciding_reads`; no "vote" left in the module.

### Mi1 — `--explain-loci` rows
- **Final status:** Deferred — needs a new explain row kind; the leftover-derived
  `reads_on_no_candidate` is now the merge's, so the remaining gap is the per-allele read counts.

### Mi2 — Locus length · **Applied**: `LocusLen::of_region(observation.region)`; test fixtures now give the locus a region as long as its reference.
### Mi3 — Flank-deletion slicing · **Applied**: `partial_read_fits_takes_the_alleles_own_flank_deletion_from_the_far_side`; mutant (right deletion ignored) killed by it alone.
### Mi4 — Error mass · **Applied**: the stutter test asserts `q_sum` −31.0 and −2.5, with the two partial rows given different masses.
### Mi5 — No crossing read · **Applied**: stated in §8.1; `a_sample_without_crossing_reads_nominates_on_two_deciding_partial_reads`.
### Mi6 — No count of what the rule did · **Deferred**: needs a run-log counter; no policy yet on where.
### Mi7 — `..AlleleSupport::default()` · **Deferred**: nothing reads those fields from the copy; a comment at the site says why only two fields are set.
### Mi8 — Older tests pass by fixture accident · **Deferred**: their claim still holds for partial reads that fit nothing; worth a sentence in each.
### Mi9 — Stale docs · **Applied**: `PartialObservation::held_flank` and `Classified::Observed::held_flank` now say the side is recorded at segments and bundles.
### Mi10 — Names · **Applied**: `far_flank`, `bases_in_allele`, `bases_in_flank`, `held_flank`, `only_fitting_allele`.
### Mi11 — Duplicated far-flank slicing · **Deferred**: sharing it means a method on `FlankDeletion` used by the likelihood too; out of scope.
### Mi12 — Long function · **Applied**: `add_deciding_reads` extracted with the order invariant in its doc.
### Mi13 — Discovery pre-pass · **Applied**: §8.1 lists it.
### Mi14 — Mate overlap · **Deferred**: Needs verification; not checked in this run.
### Mi15 — Edge-deletion exclusion · **Deferred**: open question 2. A mutant that blocked only left deletions is now killed by the right-deletion case.
### Mi16 — Public module · **Superseded** by M6 (module now private).
### Nits — **Applied**: redundant length guard removed (a boundary test shows a read past the stored flank still fits nothing); dead `sequence_may_run_on_past` guard documented; edge-deletion test covers both sides; env-reader test moot. **Not applied**: flat list instead of nested `Vec`, `super::super::fixtures` path.

## 5. Deferred findings to carry forward
- Mi1, Mi6, Mi7, Mi8, Mi11, Mi14, Mi15 — reasons above. M2 and M3 carry follow-ups (open question
  1; the segment read model's far flank).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing in `benches/` reaches candidate selection. The 50-sample, 5× run took
  107 s to walk and call, against no recorded baseline for this exact pair; not a measurement.

## 10. Commands run
All through `scripts/dev.sh` (Apple `container`), from this worktree:
- `cargo test --release --no-fail-fast -j 2 --lib allele_candidates` (before, after, and once per mutant — `tmp/review/mutants.sh`)
- `cargo fmt` / `cargo fmt --check` on the touched files
- `cargo clippy --release -j 2 --lib --tests`
- `cargo test --release --no-fail-fast -j 2 --lib run::`, `--lib calling::`, `--lib cross_platform_digests`
- `DEV_EXTRA_MOUNT=…/pop_var_caller-str-bundles/tmp/bundle_bench/runs scripts/dev.sh tmp/verify/run.sh` (build `pvc`, walk and call `pcr_1s_5x` and `pcr_50s_5x`, compare with `tmp/pc/runs/<run>/x_on.vcf`, delete the psps)
- `python3 tmp/verify/fields.py`, `tmp/verify/dpdir.py` (which fields differ, and which way `DP` moved)

## 11. Command results
- allele_candidates before fixes → 0, 188 passed
- allele_candidates after → 0, 195 passed
- mutants (5) → each fails exactly its intended test: leftover from the copy → `the_leftover_holds…`; side from the recorded flank only → `…lifts_a_once_seen_contraction`; votes on reversed samples → `deciding_reads_count_only_for_the_sample…`; right deletion ignored → `partial_read_fits_takes_the_alleles_own_flank_deletion…`; only left deletions block → `no_partial_read_counts_where_an_allele_deletes_across_the_edge`
- fmt → 0; clippy → no finding in touched files
- run:: → 0, 779 passed; calling:: → 0, 1248 passed, 3 ignored; cross_platform_digests → 0, 4 passed
- pcr_1s_5x → 1,916 records, VCF body byte-identical
- pcr_50s_5x → 7,829 records; 6,675 identical; 1,154 differ only in `INFO/DP` and `FORMAT/DP`, all lower (1,653 cells, −1,668 reads)

## 12. Notes
- Owner-required deviations from the review skill: no git worktrees for reviewers (read-only from
  this worktree), and eleven categories dispatched to four reviewers.
- HG002 was not re-run: no review fix changed a genotype, which the two benchmark runs show.
