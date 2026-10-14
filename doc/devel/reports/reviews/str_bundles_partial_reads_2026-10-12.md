# Code Review: str_bundles_partial_reads
**Date:** 2026-10-12
**Reviewer:** rust-code-review skill (orchestrator), with the rust-performance-review skill for §6.1
**Scope:** step 4b of `str_bundle_loci.md` — no periodicity check at bundles, reads that end inside a bundle scored, and the sequencing-indel term (uncommitted diff on `str-bundles`, HEAD 1ef679aa)
**Status:** Request-changes (the cost Major and the test Majors are fixed in the same session; see [fixes_applied_str_bundles_partial_reads_2026-10-12.md](fixes_applied_str_bundles_partial_reads_2026-10-12.md))

---

## 1. Scope

- **What was reviewed:** the uncommitted diff of branch `str-bundles` over 1ef679aa (9 files,
  +1,227 / −192), frozen before any fix as
  `tmp/review_2026-10-12_str-bundles-partials/step4b_under_review.patch` with full copies under
  `snapshot/`. Design: [`str_bundle_loci.md`](../../ng/spec/str_bundle_loci.md) §6.1, §7.2, §7.3,
  the new §7.6, §8, §10 (step 4b) and §11.5.
- **Reviewed against:** 1ef679aa plus the working tree as found.
- **In scope:** `select_ssr_bundle` without the periodicity check
  (`calling/allele_candidates/ssr.rs`); `bundle_censored_emission`, `CensoredLayout`, `HeldFlank`,
  the indel term in `bundle_emission` and `part_repeat_emission` (`calling/likelihood/ssr_emission.rs`);
  the trait's new `witnessed_positions` parameter and `fill_emissions` (`calling/likelihood/ssr.rs`);
  `call_one_ssr_bundle_locus` passing every observation (`run/callers.rs`);
  `BundleCandidate::reference_hull_len` (`calling/inference/repeat_tract_parameters.rs`); doc
  comments in `locus_generation/ssr.rs`, `ssr_bundle.rs`, `window_coverage/depth.rs`.
- **Out of scope:** segment scoring (unchanged; one pre-existing defect is recorded at M4 because
  the brief asked for it).
- **Categories dispatched** — read-only sub-agents in this worktree, by the owner's rule **no
  worktrees, no builds**, so no sub-agent ran a mutation; the lead ran the two mutations quoted in
  §8:
  - reliability, with the "challenge the tests" pass and a Python port of the run-by-run sum
    checked against a brute-force sum (`model.py`, `either.py`);
  - errors + idiomatic + smells + refactor safety;
  - naming + module structure + defaults + extras + float portability, and step 8a (the diff's own
    numbers).
  - Performance: by the lead, with `perf record -e cpu-clock` in the container, atomic operation
    counters and single-thread CPU times (§6.1). The skill's per-category perf fan-out was not run:
    one measured profile named the sites.
  Per-category files: `tmp/review_2026-10-12_str-bundles-partials/`.

## 2. Verdict

Request-changes. The new sum is correct — a port of it agrees with a brute-force sum over every
slip set to 4 parts in 10¹⁶ on every read shape tried — but it cost about 25 ms more per bundle at
50 samples, which does not scale to thousands; two of the new tests cannot fail; and two modelling
properties of the partial-read emission contradict or go beyond what §7.6 says (M2, M3), which are
the owner's to settle. Nothing here moves a segment.

## 3. Execution status

- `cargo fmt --check` on the touched files (rustfmt, in the container) → exit 0 after the fixes.
- `cargo clippy --release --all-targets` → exit 0; warnings on the touched files after the fixes:
  only `src/locus_generation/ssr.rs:2701` `type_complexity`, which is pre-existing.
- `cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run::callers
  window_coverage:: cli::cross_platform_digests` → exit 0, `1774 passed; 0 failed; 4 ignored`,
  including the four cross-platform digest tests and the one-run regression
  (`a_one_run_bundle_scores_as_the_segment_it_is`).
- Not run: the full `--all-targets` suite, `cargo doc`, `cargo audit` (the brief asked for touched
  tests only).
- Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. **Should a read be neutral about the runs it never reached?** (M2.) Either every unreached run
   is a factor of one, or §7.6 records that each contributes its total slip probability, which
   differs by count (measured below). Owner's decision.
2. **Should a read whose witness covers the whole hull be scored from the flank it anchored?** (M3.)
   That needs the anchored side carried from the read classifier to the emission, a change to the
   observation or the witness. Owner's decision.
3. **The segment fix** (M4) is a separate change; it moves segment calls.

## 5. Top 3 priorities

1. **H1–H3 — cost.** 25 ms more per bundle at 50 samples, linear in samples. Fixed: 8–10 ms, byte-identical calls.
2. **M5, M6 — tests that could not fail.** Fixed: both indel-term tests rewritten; three new partial-read tests; two mutations shown caught.
3. **M2, M3 — what a partial read says.** Measured; left to the owner.

## 6. Findings

### 6.1 Performance (rust-performance-review)

**Hot-path evidence.** Single thread, 50 samples, PCR, 30×, the 500 bundle hulls of
`tmp/step4b/cost/bundles.bed` ±60 bases (the implementer's set), true stutter, new psps under
`tmp/step4b_review/cost/`. The host was loaded by other work during part of the session (load
average 27), so wall times moved by a factor of two between runs; the comparisons below are CPU
times from runs made back to back in one container session.

`perf record -e cpu-clock -F 999 -g` of the build under review, self time (27K samples):

```text
26.41%  ssr_emission::substitution_probability
17.38%  libm::math::log::log                      ← all from substitution_probability
 8.37%  CensoredLayout::probability_on_from_fixed_piece
 7.47%  __memcpy_generic
 6.29%  ssr_emission::place_slip_in_each_run      ← CensoredLayout::lay_out, once per read per candidate
 5.98%  (the renderings' fold in probability_on_from_run)
 5.94%  ssr_emission::censored_from_one_flank
 4.46%  run_parse::parse_into_runs
```

Operation counts per censored call (one read, one candidate, one flank; 1.6 million calls): about
200 letter comparisons of 11.6 bases each, the first 2 of them agreeing on average; 25 run states
visited, 13 answered from the memo; 13 part-repeat places. Of 1.6 million partial emissions,
707,000 held the left flank, 708,000 the right and 185,000 the whole hull (scored twice).

**H1: `src/calling/likelihood/ssr_emission.rs` `censored_from_one_flank` — an allele's layout is rebuilt for every read** — *Major*
- **Confidence:** High (profile above; with layouts kept, `place_slip_in_each_run` fell from 6.3 to
  2.3 in 100 of a profile whose samples fell from 27K to 17K).
- **Problem:** `run_slip_ranges` and `CensoredLayout::lay_out` — every slip of every run rendered with
  `place_slip_in_each_run`, which allocates per placement, plus every part-repeat rendering — ran once
  per partial observation per candidate (twice for a whole-hull read), though none of it depends on the read.
- **Fix (applied):** keep up to 16 layouts in the row's scratch, each found by comparing everything it
  was built from (the allele's bases, the flank beyond the hull, the runs, the runs' places and counts,
  the run models, the flank held) — never an address; and score a candidate's partial reads together
  (`fill_emissions` loops candidates outside, reads inside; each cell is written once, so the order
  changes no value). Hit rate 89 in 100.

**H2: `substitution_probability` builds a `FlatEmission` per comparison — two logarithms each** — *Major*
- **Confidence:** High (17% of the profile was `libm::log`, all under this function).
- **Fix (applied):** `substitution_scores` looks the two scores up; `letters_under` compares. The
  censored path looks them up once per read. The complete-read path is unchanged.

**H3: a fixed piece reached at one offset by two paths is compared twice** — *Likely*
- **Fix (applied):** fixed-piece states are memoised beside run states.

**L1: comparisons re-add the bases that agree** — *Likely*, small
- **Fix (applied):** `MatchedPrefixLetters` keeps the sums of `i` match scores and continues from the
  first mismatch, bit for bit `letters_under`'s result (`a_matched_prefix_gives_letters_under_bit_for_bit`).
  It saves about 17 in 100 of the bases summed; the gain was within run-to-run noise.

**Result.** CPU time, single thread, the same session:

| build | calling, s (two runs) | more than step 4, per bundle |
|---|---|---|
| step 4 as reviewed (no partial reads) | 2.0, 1.9 | — |
| step 4b under review | 14.5, 13.9 | about 24 ms |
| reviewed (H1–H3, L1) | 7.0, 5.7 | about 8–10 ms |

The calls are byte-identical before and after on every run measured (§10 below). **What is left** is
the number of comparisons — about 200 per read per candidate per flank, a dependent chain of
additions each — and `parse_into_runs` (10 in 100 of the profile), which comes from candidate
selection now parsing every bundle's alleles where the periodicity check used to refuse some first.
The next lever is sharing run states between candidates whose layouts agree from a run on (the
outer run differs, the rest is the same); it was not built.

**Note:** the cost is linear in partial reads × candidates, so at a thousand samples it is about
twenty times the 50-sample figure, as step 4's is.

### 6.2 Major

**M2: `ssr_emission.rs` `CensoredLayout::lay_out` (`unreached_mass`) — a partial read is not neutral about runs it never reached** — *Major*, Deferred to the owner
- **Categories:** reliability.
- **Confidence:** High on the arithmetic (Python port; the Rust is unchanged in this respect).
- **Problem:** an unreached run contributes its total probability over the slips it can make. A run
  with one copy cannot contract, so it loses the contraction share; a run with no copy is exactly one.
  A read `GG (CA)8 G` stopping before the `TTG` run scores, under the shipped constants, 0.863 against
  no `TTG` copy, 0.811 against one, 0.852 against two and 0.854 against five or more: about 6 in 100
  between one copy and none, always in the same direction. §7.6 says such a read "says nothing" about
  those runs. The same mass is missing from a complete read's emission.
- **Fix (question):** a factor of one for every unreached run, or keep it and record the limitation.
  The spec now records the measurement (§7.6, "Found in review").

**M3: `ssr_emission.rs` `HeldFlank::of` — a read whose witness covers the whole hull is scored at half weight** — *Major*, Deferred to the owner
- **Categories:** reliability; naming and errors raised the same type at Minor.
- **Confidence:** High that the information is lost; Medium on what it moves.
- **Problem:** the classifier (`locus_generation/ssr.rs`, `partial`) knows the anchored border, but the
  witness clamps a reach of the whole hull to the whole hull, so `HeldFlank::of` returns `Either` and
  the emission averages both readings. The wrong reading is about 4 parts in 10¹³ of the right one, so
  the score is half the left-only score — up to ln 2 = 0.69 nats against the outlier floor. 185,000 of
  1.6 million partial emissions at the cost benchmark.
- **Fix (question):** carry the anchored side to the emission. The ½/½ average is a proper mixture
  for a read whose side is truly unknown; it is not double counting.

**M4 (pre-existing, segments): `ssr_emission.rs` `censored_emission` — a right-flank partial read is compared with the start of the tract** — *Major*, for a separate change
- **Confidence:** High (code read; measured).
- **Problem:** the generator stores a right-anchored partial's bases as the end of the tract
  (`region_seq[tract]`, laid from the right border), and `censored_emission` compares them with the
  first bases of the candidate. Concrete case: a `CA` tract read from the right showing `ACACA` is
  compared with `CACAC` — five mismatches, about 4 × 10⁻¹⁸ against every candidate — so the read falls
  to the outlier floor. On one simulated sample (PCR, 30×; `ng_ssr_loci_dump` and
  `tmp/step4b_review/segment_right_partials.py`), 14,175 of 22,661 right-flank partial reads at
  segments are compared out of phase: none at period 1, 1,699 of 3,412 at period 2, and 4,099 of
  5,349 at period 6 — 27 in 100 of the 51,818 partial observations. At a pure tract the loss is
  evidence thrown away rather than evidence for the wrong allele; at an interrupted tract it can be
  the second.
- **Fix:** compare a right-held read against the candidate's last bases (the bundle path's mirror
  does exactly that). It moves segment calls, so it needs its own measurement.

**M5: two indel-term tests could not fail** — *Major*, fixed
- `a_change_one_run_makes_alone_takes_no_indel_term`: a wrongly added term was 5 parts in 10⁴³ of the
  expected value. `a_change_no_single_run_makes_is_also_priced_as_a_sequencing_indel`: the slip sets
  were 1 part in 10¹⁷, so "indel term only" passed.

**M6: the partial-read tests left the memo, three read shapes and the whole-hull reading untested** — *Major*, fixed
- With two runs no memo entry is read back; reads ending in fixed sequence, at a boundary, or longer
  than any carrier were missing; the whole-hull check's right reading was below the tolerance.

**M7: nothing tests that a bundle's partial reads reach calling** — *Major*, Deferred
- The removed filter in `run/callers.rs` could return with the unit suite green; only the §10 table
  would show it. A test needs a driver-level bundle fixture.

### 6.3 Minor

- **Mi1:** `censored_from_one_flank` — the `leading_flank_bases == 0` precondition was a
  `debug_assert`; a release build would misplace the runs silently. Fixed: `assert_eq!`.
- **Mi2:** `HeldFlank::Either` reaches helpers that mean Left or Right (`flank_beyond_the_hull` and the
  part-repeat comparison read it as Right, `lay_out` as Left), and a witness touching neither edge
  falls into `Either`. Not reachable today. Deferred (with M3).
- **Mi3:** `BundleCandidate::reference_hull_len` is a bare `u32` set with `as u32` where `LocusLen`
  exists, and nothing checks it equals the witness's locus length. Deferred.
- **Mi4:** `flank_beyond_the_hull` returns an empty flank when the deletion is longer than the flank,
  where the complete path's `respelled_over` returns `None`. Deferred.
- **Mi5:** the `unreachable!` in `fill_emissions` could be removed by having
  `partial_observations()` yield the positions. Deferred.
- **Mi6:** stale documents — `locus_generation/ssr.rs` `MeasuredRepeat::bundle_hull` ("calling leaves
  it out"), `locus_generation_ssr.md` line 116, `str_bundle_loci.md` §6.1 (`call_one_ssr_bundle_locus`
  "complete observations only") and §7.2's three periodicity bullets, and `max_off_grid_share`'s doc.
  Fixed.
- **Mi7:** `lay_out` is long and does three jobs; `widest = 2·MAX_PART_REPEAT_SLIP` restates the
  model's bound by hand (correct: the widest priced change is 19 bases, at period 2); the scratch
  buffer tuples. `lay_out`'s buffers are now one named struct; the rest Deferred.

### 6.4 Nits

"About 0.0023" in §7.3 is the (−1, +1) slip set alone; the full sum is 0.00237 (now 0.0024). A test
doc said reads start "from two bases"; for `AAG` they start at three (fixed). The step 4b table's
generic-path row covers 10,457 variant cells against the column's 10,435 (that caller's scoring
assigns cells differently; not changed). `ssr_emission.rs` is now about 4,900 lines.

## 7. Out of scope observations

- M4 above, filed as a finding because the brief asked for it.
- `locus_generation/ssr.rs:2701` `type_complexity` (clippy), pre-existing.

## 8. Missing tests to add now

Added in this session (all pass):

- `a_partial_read_ending_in_fixed_sequence_or_at_a_boundary_sums_every_rendering` — reads ending
  inside `GT`, at the `CA` run's end, at `GT`'s end, from each flank, against (8,5), (9,1), (3,6);
  and a read longer than any carrier, which scores exactly zero.
- `a_partial_read_over_three_runs_sums_every_slip_set` — `G (CA)6 TT (AAG)4 CC (TG)5 A`, from each
  flank. **Mutation:** the run memo keyed by offset alone fails it, and passes every earlier test.
- `an_either_read_averages_both_readings_where_both_explain_it` — `(CA)10` between flanks that
  continue the run. **Mutation:** dropping the right reading fails it.
- `a_matched_prefix_gives_letters_under_bit_for_bit` — 3,444 comparisons, bits equal.
- `a_change_no_single_run_makes_is_also_priced_as_a_sequencing_indel` and
  `a_change_one_run_makes_alone_takes_no_indel_term`, rewritten so each half is pinned by a factor of
  a million.

Still missing: `a_bundle_sample_with_only_partial_reads_gets_a_non_flat_likelihood` (M7).

## 9. What's good

- The run-by-run sum is right on every read shape tried, including the memo and both flanks
  (`ssr_emission.rs`, `CensoredLayout`).
- `part_repeat_emission` takes the letters as a closure, so the complete and partial reads price
  the same renderings (`ssr_emission.rs`); the fix built the partial reads' renderings from the same
  place walker.
- The one-run regression against the segment's own censored sum, over pure, interrupted and
  mononucleotide tracts (`a_one_run_bundles_partial_read_scores_as_the_segments`).

## 10. Commands to re-verify

- `scripts/dev.sh bash tmp/step4b_review/check.sh` — fmt check and clippy on the touched files, the
  touched modules' tests.
- `scripts/dev.sh bash tmp/step4b_review/cost.sh <label> <binary> 2` — the cost benchmark.
- `DEV_EXTRA_MOUNT=<main>/benchmarks scripts/dev.sh bash tmp/step4b_review/bench.sh` — new psps,
  pre-review and reviewed calls on the three benchmark runs, VCF bodies compared, cells scored;
  HG002 30× compared with `tmp/step4b/hg002/step4b_30x.vcf` and scored per region.

## 11. Measurements after the fixes

**The calls did not change.** VCF bodies (everything but the `##` header lines) of the reviewed
binary against the binary under review, new psps for each run:

| run | records | result |
|---|---|---|
| 1 sample, PCR, 5× | 1,898 | byte-identical |
| 1 sample, PCR, 30× | 2,591 | byte-identical |
| 50 samples, PCR, 30× | 6,005 | byte-identical |
| HG002 30× (call-from-alignments, against `tmp/step4b/hg002/step4b_30x.vcf`) | 20,347 | byte-identical |
| cost benchmark (500 bundles, 50 samples) | — | byte-identical |

Bundle cells with a variant, right / wrong / missed; false calls — as the implementer reported:
1 sample 5× 99 / 29 / 80, 0; 1 sample 30× 207 / 1 / 0, 0; 50 samples 30× 10,400 / 13 / 22, 1;
HG002 30× per region 862 / 150 / 119, 20. Segment cells against main: 0 of 3,000, 0 of 3,000 and
0 of 150,000 differ.

**Cost**, single thread, CPU seconds of calling, two sessions under different background load:

| build | session A | session B |
|---|---|---|
| step 4 as reviewed (no partial reads) | 2.0, 1.9 | 3.0, 2.2 |
| step 4b under review | 14.5, 13.9 | 23.8, 14.3 |
| reviewed | 7.0, 5.7 | 9.4, 6.6 |

Partial reads add about 8–9 ms a bundle at 50 samples where they added about 24 ms — about three
times less; calling all 500 bundles still takes about three times step 4's time.
