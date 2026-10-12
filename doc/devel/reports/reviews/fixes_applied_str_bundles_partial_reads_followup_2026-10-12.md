# Fix Application Report: str_bundles_partial_reads_followup_2026-10-12.md

**Date:** 2026-10-12
**Source review:** `doc/devel/reports/reviews/str_bundles_partial_reads_followup_2026-10-12.md`
**Source state reviewed against:** `str-bundles` at 430eb4ed plus the uncommitted step 4c diff
**Execution mode:** non-interactive (session lead's brief; nothing committed)
**Overall status:** Completed — every Blocker and Major applied; Minors partly applied, the rest deferred

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 4
- Minors: 13
- Nits: 6 groups

### Outcome totals
- Applied: 5 (B1, M2, M3, Mi1, Mi13)
- Applied with adaptation: 2 (M1, M4)
- Already fixed: 0
- Deferred: 11 (Mi2 — its test half applied — and Mi3–Mi12); the Nits too
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### What the two required fixes did

- **The substitution rate (M1).** A partial read at a bundle now compares each run's bases at that
  run's own rate — the rate fitted for the run's count in that candidate — and every other base
  (fixed pieces, the flank beyond the hull) at one rate per read group that no candidate moves: the
  mean of the runs' rates at their *reference* counts. A complete read keeps the cell's mean rate.
  - Proof of neutrality: on `GG (CA)8 <spacer> (TTG)n C` with rates of 10⁻³·(1 + count), a read
    stopping in the spacer scores bit for bit alike against 1, 2, 5, 9 and 12 `TTG` copies, from
    either flank; at the old mean rate the same read scores more than 5% apart against 5 and 12.
  - The one-run regression is unchanged: one run's rate is the cell's rate exactly, and
    `a_one_run_bundles_partial_read_scores_as_the_segments` passes at its twelve digits.
  - Effect on the benchmarks: **none**, as expected — every rate on those runs is the defaulted
    constant, and the VCFs were byte-identical before and after (§9). It moves calls only on a run
    with fitted rates.
- **The cost (M4).** The partial-read walk now goes forward, carrying the probability of reaching
  each state (a piece entered with so many of the read's bases spent), and at bundles of two runs or
  more drops states within a total budget of 10⁻¹² per read (`NEGLIGIBLE_PARTIAL_READ_MASS`): the
  read comes out at most 10⁻¹² low, never high. One-thread CPU on the 500-bundle cost set, three
  alternating runs: **6.8–6.9 s, against 16.4–17.1 s before and 12.0–12.3 s for step 4b** — and the
  VCFs byte-identical to before on all three benchmark runs, HG002 30× and the cost set.

### Validation summary
- `rustfmt --edition 2024 --check` on the 38 touched `.rs` files → exit 0, clean.
- `cargo clippy --release --all-targets` → exit 0; nothing on a line this run touched. The warnings
  left are earlier ones (`examples/dhat_ng_*`, `src/run/edge_deletion_gathering.rs:201`).
- `cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run:: window_coverage:: psp:: cli:: parameter_estimation::`
  → exit 0, `3637 passed; 0 failed; 11 ignored`, the four `cli::cross_platform_digests` tests
  included.
- `cargo test --release --no-fail-fast --test ng_calling_loop_calls_genotypes` → exit 0, `16 passed`.
- `cargo test --release --lib -- ssr_emission:: psp::record::` after the last edit → exit 0,
  `149 passed`.
- `cargo doc --no-deps`, `cargo audit` → not run: no public API or dependency change beyond the new
  `pub` constant and context fields, whose doc links clippy's doc lints accepted.
- Performance check → run as the review's own measurement (§9), not as a criterion comparison: no
  bench in `benches/` reaches the bundle emission.

### Unresolved high-priority findings
- None. Mi2 (an older build may extend a newer minor's psp) is the most consequential deferral;
  its only path is `PspWriter::append`, which has test-only callers.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | No test fails if a right-held read records the wrong flank | Apply | Applied | No | `src/locus_generation/ssr.rs` | Pass; mutation fails it | No |
| M1 | Major | Partial read's letters depend on an unreached run through the rate | Apply | Applied with adaptation | No | `ssr_emission.rs`, `repeat_tract_parameters.rs`, spec | Pass | No |
| M2 | Major | Row's hand-off of `held_flank` untested | Apply | Applied | No | `src/calling/likelihood/ssr.rs` | Pass | No |
| M3 | Major | "Price the same renderings" false; partial ≥ complete unpinned | Apply | Applied | No | spec, `ssr_emission.rs` | Pass | No |
| M4 | Major | Step costs 1.37× step 4b | Apply | Applied with adaptation | No | `ssr_emission.rs`, spec | Pass; VCFs identical | No |
| Mi1 | Minor | Decoder accepts duplicates and a flank on a complete observation | Apply | Applied | No | `src/psp/record.rs` | Pass | No |
| Mi2 | Minor | `check_this_build_writes` accepts extra trailing fields | Defer | Deferred (test half applied) | No | `src/psp/record.rs` (test) | Pass | Yes |
| Mi3 | Minor | `held_flank` TSV column untested | Defer | Deferred | No | None | N/A | Yes |
| Mi4 | Minor | Partials' sort key written twice | Defer | Deferred | No | None | N/A | Yes |
| Mi5 | Minor | Format version as four bools | Defer | Deferred | No | None | N/A | Yes |
| Mi6 | Minor | `CensoredSlip` names indels too | Defer | Deferred | No | None | N/A | Yes |
| Mi7 | Minor | Indel-place rule written twice | Defer | Deferred | No | None | N/A | Yes |
| Mi8 | Minor | Hand-counted TSV padding | Defer | Deferred | No | None | N/A | Yes |
| Mi9 | Minor | "The review's M2/M3" names no review | Defer | Deferred | No | None | N/A | Yes |
| Mi10 | Minor | `held_flank` representable on a complete observation | Defer | Deferred (decoder refuses it, Mi1) | No | None | N/A | Yes |
| Mi11 | Minor | `ObservationIdentity` tuple | Defer | Deferred | No | None | N/A | Yes |
| Mi12 | Minor | Silent ½/½ fallback | Defer | Deferred | No | None | N/A | Yes |
| Mi13 | Minor | HG002 column not over shared psps | Apply | Applied | No | spec §10 | N/A (doc) | No |
| Nits | Nit | six groups | Defer | Deferred | No | None | N/A | No |

## 3. Questions asked and answers

None asked. The two open questions of the review were settled by the session lead's brief: the rate
rule was to be "the least invasive correct rule", and cost levers "pruning negligible indel branches
with a stated bound" were named. The choices made are in the per-finding log and are the owner's to
overrule.

## 4. Per-finding log

### B1 — No test fails if a right-held read records the wrong flank
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** the code is right; the gap is the test, and a wrong flank scores 25,214 reads at
  50 samples 30× from the wrong end.
- **Implementation summary:** `a_partial_over_the_whole_hull_from_the_right_records_the_right_flank`:
  five extra `CA` bases where the left flank would be, the read placed three bases into the left
  flank, anchoring only the right; it must be a partial over the whole hull recording
  `Some(HeldFlank::Right)`.
- **Review suggestion used verbatim?:** No (the reviewer named the test; the fixture is new).
- **Adaptation:** None.
- **Verification performed:** mutation `.then_some(border)` → `.then_some(HeldFlank::Left)`: the new
  test fails (`left: Some(Left)`, `right: Some(Right)`), the left-held test still passes; reverted.
- **Files changed:** `src/locus_generation/ssr.rs`
- **Tests added or modified:** the one above.
- **Validation:** `cargo test --release --lib -- whole_hull` → exit 0, 2 passed (101 under the mutation).
- **User input:** None. **Follow-up:** None. **Residual risk:** None.

### M1 — A partial read's letters depend on an unreached run through the substitution rate
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** the same defect M2 fixed for stutter, by another route.
- **Implementation summary:**
  - `SsrScoringContext` gains `run_substitution_rates` (one per run, hull order) and
    `fixed_substitution_rate`; `SsrScoringContext::with_partial_read_rates` sets them on a bundle
    context. A context built without them (segments, the tests that predate them) compares a partial
    read at the cell's rate throughout, as before.
  - `gather_for_bundle` records each run's rate as it already looked it up (a run with no copy has
    no bases and takes the defaulted constant unread), and one fixed rate per read group: the mean of
    the runs' rates at their reference counts. That rate's warrant and its falling back are folded
    into every cell of the read group, as the cell's own rates are.
  - The censored walk keeps one `MatchedPrefixLetters` per distinct rate (the fixed one first, then
    each run whose rate differs), each keyed on its rate's bits, so on a run with one rate everywhere
    it is one table as before.
- **Review suggestion used verbatim?:** N/A (the brief named the shape: per-run rates on the runs'
  bases, fixed sequence at a defined rate).
- **Adaptation:** the "defined rate" for fixed bases is the reference counts' mean, per read group —
  count-independent, so neutral, and still the read group's own fit. A complete read is unchanged.
- **Verification performed:** the neutrality test below; its second half reproduces the defect at
  the mean rate (more than 5% apart against 5 and 12 copies). The one-run regression passes at
  twelve digits. VCFs byte-identical on every benchmark run (all rates defaulted there).
- **Files changed:** `src/calling/likelihood/ssr_emission.rs`,
  `src/calling/inference/repeat_tract_parameters.rs`, `doc/devel/ng/spec/str_bundle_loci.md` §7.6.
- **Tests added:** `a_partial_read_is_compared_at_the_rates_of_what_it_reached`,
  `a_bundle_cell_carries_each_runs_rate_and_a_fixed_rate_no_candidate_moves`.
- **Validation:** `cargo test --release --lib -- ssr_emission:: repeat_tract_parameters:: likelihood::ssr:: summarise_condition::` → exit 0, 236 passed.
- **User input:** None. **Follow-up:** None.
- **Residual risk:** a partial and a complete read of the same bases are now compared at slightly
  different rates where rates are fitted per count. Both are honest; they differ by the spread of
  the runs' fitted rates.

### M2 — The row's hand-off of `held_flank` is untested
- **Severity:** Major
- **Final status:** Applied
- **Implementation summary:** `the_row_hands_the_model_each_partials_held_flank`: a model that
  records the flank of every censored call, wrapped around the real one; two partials of the same
  bases and witness, one recorded `Right` and one with none, must reach it as exactly those, once
  per candidate.
- **Files changed:** `src/calling/likelihood/ssr.rs`
- **Validation:** passes; passing `None` at the call site would give four `None`s and fail it.
- **Follow-up / residual risk:** None.

### M3 — "The two rules price the same renderings" is false, and nothing pins how they relate
- **Severity:** Major
- **Final status:** Applied
- **Implementation summary:** §7.6 corrected: what differs (a partial lets every reached run take an
  indel and adds them; a complete read takes one only off every single run's slip lengths, averaged)
  and what holds (each run's changes a distribution with no length counted twice; partial ≥
  complete). New test `a_two_run_partial_read_is_at_least_as_likely_as_the_complete_read`: 18 reads
  — the reference hull, one run a unit longer or shorter in four ways, each with a base inserted or
  removed — scored complete and as a left-held partial over the whole hull.
- **Files changed:** `doc/devel/ng/spec/str_bundle_loci.md`, `src/calling/likelihood/ssr_emission.rs`
- **Validation:** passes. **Follow-up:** whether the complete read should adopt the per-run rule is
  a modelling question for the owner, not taken here.

### M4 — The step costs 1.37 times step 4b
- **Severity:** Major
- **Final status:** Applied with adaptation
- **Profile before** (`perf record -e cpu-clock`, profiling build, one thread, 500 bundles, 50
  samples PCR 30×): censored emission 72% inclusive; `MatchedPrefixLetters::letters` 28% self, the
  averaging closure around it 18%, `probability_on_from_run` 9.5%, `parse_into_runs` 8.3%.
- **Implementation summary:** `probability_on_from_fixed_piece` and `probability_on_from_run` (a
  recursive, memoised walk) replaced by `walk_from_the_held_flank`, which goes forward: states in
  piece order, offsets ascending, each carrying the summed probability of every path reaching it.
  From two runs on, a state no more probable than what is left of a 10⁻¹² budget is dropped. Every
  factor after a state is at most one, so the read comes out at most 10⁻¹² low and never high; the
  hand-sum tests now check exactly that (`matches_the_hand_sum`: not above beyond rounding, not
  below by more than the budget). A one-run bundle is walked whole.
- **Adaptation:** of the three levers the brief named, only a bounded prune was used, and it prunes
  states rather than indel branches, since a state is where the cost sits. The two others were not
  needed.
- **Measured** (one thread, the cost set, three alternating runs each, user CPU):

  | build | CPU, s |
  |---|---|
  | step 4b, reviewed | 11.95–12.30 |
  | step 4c as reviewed | 16.37–17.11 |
  | step 4c with the rate fix only | 16.29–18.73 |
  | forward walk, nothing dropped (experiment, two runs) | 14.06–14.41 |
  | step 4c reviewed (both fixes) | 6.76–6.91 |

  So 0.41 times the step 4c cost and 0.57 times step 4b's — about 20 ms less a bundle locus. The
  rows were measured in two sessions under the same load; the ratios within a session are the
  measurement.
- **Byte identity:** the VCF and the parameters file of the cost set are identical across step 4c,
  the rate fix, the forward walk with and without the budget, and the final build (headers naming the
  command line excepted). VCF bodies identical between step 4c and the final build on pcr_1s_5x,
  pcr_1s_30x, pcr_50s_30x and HG002 30×.
- **Files changed:** `src/calling/likelihood/ssr_emission.rs`, `doc/devel/ng/spec/str_bundle_loci.md`
  §7.6 and §10.
- **Residual risk:** the budget is an absolute 10⁻¹² per read. A read whose every allele explains it
  worse than that can come out as zero instead of a smaller number; the row's outlier floor near
  0.2 / 40 takes such a read either way.

### Mi1 — Decoder leniency
- **Final status:** Applied. A repeated observation index and a flank on a complete observation are
  refused as malformed. The round-trip test now records flanks on two partial observations and
  damages the flank, the index (duplicate) and the index (a complete observation), each refused.
- **Files changed:** `src/psp/record.rs`. **Validation:** `psp::record::` tests pass.

### Mi2 — `check_this_build_writes` accepts extra trailing fields
- **Final status:** Deferred, test half applied: `a_format_2_0_manifest_is_read_as_carrying_no_held_flanks`
  now asserts this build refuses to extend a 2.0 file. Refusing a manifest *longer* than this
  build's needs a new `RecordLayoutError` variant, and the shipped 2.0 binary cannot be changed
  anyway; its only path is `PspWriter::append`, test-only today.

### Mi3–Mi12, Nits — Deferred
- Readability and structure (Mi4–Mi9, Mi11, Nits) or a type change touching every observation site
  (Mi10), none changing a result; Mi3 and Mi12 add a test and a run-report counter, outside the
  brief's two required items. Carried forward below.

### Mi13 — HG002 column
- **Final status:** Applied in §10 step 4c (not over shared psps; one record more than step 4b).

## 5. Deferred findings to carry forward
- Mi2 — refuse a manifest longer than this build's in `check_this_build_writes`.
- Mi3 — test the `held_flank` column's values.
- Mi4 — one `PartialObservation::order_key()`.
- Mi5 — `BodyFormat` enum for the record layout.
- Mi6 — rename `CensoredSlip` / `slips` / `slips_of_run` to `RunChange` / `changes` / `changes_of_run`.
- Mi7 — one helper for where a run can take a sequencing indel.
- Mi8 — derive the TSV padding from `COLUMNS.len()`.
- Mi9 — name the review behind "finding M2/M3".
- Mi10 — move `held_flank` into `ReadWitness::Partial`.
- Mi11 — `ObservationIdentity` as a struct.
- Mi12 — count whole-hull partials scored ½/½ in the run report.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** yes — M4 is a cost change on the calling path.
- **Baseline saved:** yes, as a binary, before any fix (`tmp/bins/pvc_step4c_prereview`).
- **Benches run:** the review's cost set (`tmp/step4c_review/cost.sh`), not `benches/`: no criterion
  bench reaches the bundle emission.
- **Outcome:** improved, 16.4–17.1 s → 6.8–6.9 s, results byte-identical.
- **Measurements, bundle cells with a variant (right / wrong / missed; false calls):**

  | | 1 sample, PCR, 5× | 1 sample, PCR, 30× | 50 samples, PCR, 30× | HG002 30× |
  |---|---|---|---|---|
  | step 4c as reviewed | 99 / 29 / 80; 0 | 207 / 1 / 0; 0 | 10,407 / 11 / 17; 1 | 862 / 151 / 118; 20 |
  | step 4c reviewed | 99 / 29 / 80; 0 | 207 / 1 / 0; 0 | 10,407 / 11 / 17; 1 | 862 / 151 / 118; 20 |

  Every STR segment cell unchanged against main (3,000, 3,000 and 150,000). Whole-hull partials
  carrying their flank, none falling back: 184 / 1,030 / 50,776, as before.

## 10. Commands run
- `scripts/dev.sh bash tmp/step4c_review/build.sh <name> [profiling]`
- `scripts/dev.sh bash tmp/step4c_review/profile.sh pvc_step4c_prereview_prof base`
- `scripts/dev.sh bash tmp/step4c_review/test1.sh`, `test2.sh` (targeted tests; the B1 mutation)
- `scripts/dev.sh bash tmp/step4c_review/cost.sh 3 …` (three sessions of alternating timings)
- `scripts/dev.sh bash tmp/step4c_review/check.sh`, `check2.sh` (fmt, clippy, tests, final build)
- `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks scripts/dev.sh bash tmp/step4c_review/bench.sh`
- `scripts/dev.sh bash tmp/step4c_review/cleanup.sh` (psps and perf data deleted; 0 psps left)

## 11. Command results
- Results as quoted in §1, §4 and §9; logs under `tmp/step4c_review/`.

## 12. Notes
- The final binary is `tmp/bins/pvc_step4c_reviewed`.
- A first run of `check.sh` hung in its container: `git diff` cannot run inside the container (the
  worktree's `.git` points at the main checkout, which is not mounted), so `rustfmt` got no files
  and waited on standard input. The script now reads the touched-file list from
  `tmp/step4c_review/touched.txt`, written on the host. That container
  (`cbd7a46b-6c7c-4037-ac31-2230ddc05f5a`, started 11:49) was still running at the end of the session;
  stopping it was not permitted from this session and is left to the owner.
