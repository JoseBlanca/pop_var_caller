# Code Review: str_bundles_partial_reads_followup
**Date:** 2026-10-12
**Reviewer:** rust-code-review skill (orchestrator), with the cost part run as a rust-performance-review
**Scope:** step 4c of the STR-bundles feature, uncommitted on `str-bundles` over 430eb4ed — the fixes for findings M2 and M3 of [str_bundles_partial_reads_2026-10-12.md](str_bundles_partial_reads_2026-10-12.md)
**Status:** Request-changes (all Blocker and Major findings were then fixed the same day; see [fixes_applied_str_bundles_partial_reads_followup_2026-10-12.md](fixes_applied_str_bundles_partial_reads_followup_2026-10-12.md))

---

## 1. Scope

- **What was reviewed:** the working-tree diff on 430eb4ed (40 files, +1,047 / −432), frozen for the
  reviewers at `tmp/review_2026-10-12_str-bundles-followup/reviewed.diff` with every changed file copied
  under `snapshot/`. Line numbers below are the snapshot's unless a finding says otherwise.
- **What the change does.** Two things, both at STR bundles (loci of several repeat runs in one hull),
  both about *partial reads* — reads that held one flank and ran out inside the hull:
  - **M2:** a partial read is scored by a walk from the flank it held, piece by piece; once its bases
    are spent the rest contributes exactly one, so a run it never reached says nothing. Each run it
    reaches offers its own slips and its own sequencing indels (length changes that are not whole
    units).
  - **M3:** the classifier records which flank a partial read held where its witness covers the whole
    hull and so cannot say it (`SequenceObservation::held_flank`). The flank is part of the tally key,
    is carried through the cohort merge and evidence shaping into the emission, and is stored in the
    psp as format 2.1's new last body field `observation-held-flanks`.
- **In-scope files:** `src/calling/likelihood/ssr_emission.rs`, `src/calling/likelihood/ssr.rs`,
  `src/locus_generation/{mod.rs,ssr.rs}`, `src/psp/{record.rs,header.rs,block.rs,walk.rs}`,
  `src/run/cohort_merge/build.rs`, `src/calling/evidence_shaping.rs`, `src/cli/inspect_psp.rs`,
  `doc/devel/ng/spec/str_bundle_loci.md` §7.6 and §10 step 4c, `doc/devel/ng/spec/psp_file_format.md`
  "Format 2.1", and the fixture-only edits in the other 25 files. Also, by the owner's request,
  `src/calling/inference/repeat_tract_parameters.rs` (`gather_for_bundle`), where the substitution
  rate a bundle cell is scored at is built.
- **Deliberately out of scope:** the complete-read bundle emission (`bundle_emission`), reviewed in
  [str_bundles_calling_2026-10-11.md](str_bundles_calling_2026-10-11.md) and unchanged here.
- **Categories dispatched** (owner-required deviation: no git worktrees; every reviewer read-only from
  the frozen snapshot, no builds, no mutation runs — the orchestrator ran the mutations it quotes):
  - `reliability` — always; with the owner's questions on psp compatibility, the tally key and the
    two indel rules.
  - `errors` + `defaults` + `extras` — always / a parser of untrusted bytes, a stable TSV, a PR; plus
    the check of every number in the diff's prose.
  - `naming` + `idiomatic` + `smells` — always.
  - `refactor_safety` + `module_structure` — always / the diff spans 40 files.
  - `float_portability` — the emission's arithmetic reaches every bundle genotype.
  - Not dispatched: `unsafe_concurrency` (no `unsafe`, no shared state added), `tooling` (no
    `Cargo.toml` change).
  - The cost question was a performance review run by the orchestrator with a sampling profile (§6, M4).

## 2. Verdict

**Request-changes.** The two fixes do what the spec says, and the psp compatibility claims hold. But
one input that decides 50,776 reads' flank at 50 samples 30× had no test that could fail (B1); the
substitution rate still carried a run the read never reached (M1); and the step cost 1.37 times the
CPU of step 4b (M4). All three, and the other two Majors, were fixed the same day.

## 3. Execution status

- Reviewers ran nothing (owner's rule). The orchestrator ran, in the dev container:
  - a `perf record -e cpu-clock` profile of one-thread calling over the 500-bundle cost set (§6, M4);
  - the timing runs, the mutation quoted under B1, and the test and benchmark runs listed in the
    fix report.
- The implementer's logs were read as they stand: `tmp/step4c/fmt.log` empty, `tmp/step4c/clippy.log`,
  `tmp/step4c/test_touched.log`.
- Findings labelled "Needs verification": none.

## 4. Open questions and assumptions

1. **Which rate should a partial read's fixed bases be compared at?** (M1.) It must not depend on any
   candidate's counts, or a read ending in fixed sequence would hear about runs it never reached. The
   fix takes the mean of the runs' rates *at their reference counts*, per read group. Any
   count-independent rate would do; this one keeps the read group's own fit.
2. **Is a dropped-state budget acceptable in place of the exact sum?** (M4.) The precedent is
   `NEGLIGIBLE_SLIP_SET_PROBABILITY` (10⁻¹²) for complete reads at three runs or more. The fix uses
   the same 10⁻¹² as a total budget per read, and only from two runs on.

## 5. Top 3 priorities

1. **B1** — a test that fails if the classifier records the wrong flank for a right-held read.
2. **M1** — score a partial read's letters at rates that depend only on what it reached.
3. **M4** — bring the step's cost down without changing a call.

## 6. Findings

### Blocker

#### B1: [src/locus_generation/ssr.rs](../../../../src/locus_generation/ssr.rs#L1753) — No test can fail if the classifier records the wrong flank for a read that held the right flank
- **Confidence:** High
- **Categories:** reliability
- **Problem:** `.then_some(border)` decides the flank of every whole-hull partial. Its one test,
  `a_partial_over_the_whole_hull_records_the_flank_it_held`, uses a left-held read only, so
  `.then_some(HeldFlank::Left)` passes every test.
- **Why it matters:** such a read would be scored from the wrong end, where it fits about 4 parts in
  10¹³ as well (the comment in `a_partial_read_sums_every_rendering_from_the_flank_it_held`), and fall
  to the outlier floor for every allele — worse than the ½/½ average this step replaced. Half of the
  50,776 whole-hull partials at 50 samples 30× held the right flank (25,214).
- **Suggested fix:** the mirror test: a read starting with five extra `CA` bases where the left flank
  would be, anchoring the right flank only, must record `Some(HeldFlank::Right)`.

### Major

#### M1: [src/calling/inference/repeat_tract_parameters.rs](../../../../src/calling/inference/repeat_tract_parameters.rs#L476) — A partial read's letters still depend on a run it never reached, through the substitution rate
- **Confidence:** High
- **Categories:** owner's question (raised by the implementer); float_portability set the constraints
- **Problem:** `gather_for_bundle` gives each cell one substitution rate, the mean of its runs' fitted
  rates, each looked up by the run's count in that candidate. The censored walk compares every base
  at that one rate. Two alleles differing only in a run the read never reached have different means,
  so the read's letters differ. Measured by the new test on `GG (CA)8 <20-base spacer> (TTG)n C`
  with a rate of 10⁻³·(1 + count): the read stopping in the spacer scores 1 in 20 or more apart
  (|ln ratio| > 0.05) against 5 and 12 `TTG` copies at the mean rate.
- **Why it matters:** it is the same defect M2 fixed for stutter, by another route; it is invisible on
  the benchmarks, whose rates are all the defaulted constant, and real on any run with fitted rates.
- **Suggested fix:** compare each run's bases at that run's rate and every other base at one rate no
  candidate moves (open question 1). A complete read keeps the cell's mean.

#### M2: [src/calling/likelihood/ssr.rs](../../../../src/calling/likelihood/ssr.rs#L876) — The row's hand-off of `held_flank` to the emission is untested
- **Confidence:** High
- **Categories:** reliability
- **Problem:** `observation.held_flank` is passed at one call site. Passing `None` there passes every
  test, and silently brings back the ½/½ average for every whole-hull partial — up to ln 2 = 0.69
  nats a read.
- **Suggested fix:** a recording model wrapped around the real one, asserting the row hands it each
  partial's own flank.

#### M3: [doc/devel/ng/spec/str_bundle_loci.md](../../../../doc/devel/ng/spec/str_bundle_loci.md#L605) — "The two rules price the same renderings" is false at two runs or more, and nothing pins how the two relate
- **Confidence:** High
- **Categories:** reliability; errors/extras (cross-category)
- **Problem:** a partial read lets every run it reaches take a sequencing indel whatever the other
  runs do, and adds them up; a complete read takes one only where no single run's slip reaches the
  length, averages it over the runs and has the others not slip. So a partial read prices renderings
  the complete rule gives zero — at `(CA)8 GT (TTG)5`, the `CA` run three bases longer, a length the
  `TTG` run's slip reaches alone.
  - What does hold, checked term by term by the reliability reviewer: each run's changes sum to at
    most one; the indels skip every whole-unit length, so no length is counted twice; and the
    partial sum holds every term of the complete one, so a partial read is never less likely than
    the complete read of the same bases. Segments pin that last property in a test; bundles did not.
- **Suggested fix:** correct the sentence, say what does hold, and add the bundle test.

#### M4: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L2053-L2138) — The step costs 1.37 times the CPU of step 4b
- **Confidence:** High (profiled)
- **Categories:** performance (owner's requirement); errors/extras and naming reviewers flagged the
  same site
- **Problem:** one-thread calling over the 500 bundle loci, 50 samples PCR 30×, on the same psps:
  16.3–16.4 s of user CPU for this build against 11.8–12.0 s for step 4b, three alternating runs each.
  The profile: the censored emission is 72% of all samples, of which comparing the read with the
  renderings is 46% (`MatchedPrefixLetters::letters` 28%, the averaging closure around it 18%) and
  the recursion itself 13%.
  - The cause: each sequencing indel shifts the read against every later piece by a non-whole-unit
    number of bases. The recursive walk then computes every state it can reach — a piece entered
    with so many of the read's bases spent — including the many where the read is already a few
    mismatches out of phase, each mismatch a factor near ε/3 ≈ 3 in 10,000.
- **Suggested fix:** walk forward, carrying the probability of reaching each state, and drop states
  too improbable to matter within a stated total budget (open question 2). The experiment that
  decides it: alternating CPU timings and a byte comparison of the VCFs.

### Minor

#### Mi1: [src/psp/record.rs](../../../../src/psp/record.rs#L1435) — The held-flank decoder accepts a duplicate entry and a flank on a complete observation
- **Confidence:** High
- **Categories:** errors/extras, reliability, refactor_safety, naming (convergent)
- **Problem:** the decoder refuses only an index past the observations and a flank other than 0 or 1.
  A repeated index silently keeps the last value. A flank on an observation whose read spanned the
  locus is accepted, though no writer can produce it. The test checks only the flank refusal, through
  a Debug string.
- **Suggested fix:** refuse both cases as malformed, and damage the index as well as the flank in the
  test.

#### Mi2: [src/psp/record.rs](../../../../src/psp/record.rs#L630) — `check_this_build_writes` accepts a manifest with extra trailing fields
- **Confidence:** High on the code path; Low on reachability (append has test-only callers)
- **Categories:** reliability
- **Problem:** a 2.0 build may append to a 2.1 file, writing records without the held-flanks field
  under a header declaring it; every reader would then refuse them as truncated. The 2.1 build has
  the same hole for a future 2.2 file. The new 2.0 test also does not pin that this build refuses to
  extend a 2.0 file, as its 1.2 sibling does.
- **Suggested fix:** require the manifest's field count to equal this build's; assert the refusal in
  the 2.0 test.

#### Mi3: [src/cli/inspect_psp.rs](../../../../src/cli/inspect_psp.rs#L317) — Nothing tests the values of the new `held_flank` column
- **Confidence:** High
- **Categories:** reliability, errors/extras
- **Problem:** swapping `left` and `right`, or printing nothing, passes. The benchmark's evidence that
  no whole-hull read fell back to the ½/½ average is read through this column.
- **Suggested fix:** dump a record with one observation of each flank and assert the last field.

#### Mi4: [src/run/records.rs](../../../../src/run/records.rs#L349) — The partials' sort key is written twice, and the copies now disagree
- **Confidence:** High
- **Categories:** refactor_safety, reliability (nit), errors/extras (nit)
- **Problem:** the merge sorts by `(witness, read group, bases, held flank)`; the re-read path by
  `(witness, read group, bases, sequence_may_run_on_past)`, under a comment saying it is the merge's
  order. Harmless today: the re-read path runs at ordinary loci only, where the flank is always
  `None`.
- **Suggested fix:** one `PartialObservation::order_key()` from an exhaustive destructure, used at both.

#### Mi5: [src/psp/record.rs](../../../../src/psp/record.rs#L568) — The format version is held as four positional bools
- **Confidence:** High
- **Categories:** refactor_safety, naming/idiomatic (convergent)
- **Problem:** four bools take 5 of their 16 combinations, one per format 1.0–2.1; swapping two
  literals compiles. The comment on the order of attempts does not match the code.
- **Suggested fix:** an ordered `BodyFormat` enum, newest first, with the bools as methods.

#### Mi6: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L631) — `CensoredSlip`, `slips`, `slips_of_run` now hold sequencing indels too
- **Confidence:** High
- **Categories:** naming
- **Suggested fix:** `RunChange`, `changes`, `changes_of_run`.

#### Mi7: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L2014) — Where a run can take a sequencing indel is decided in two places
- **Confidence:** High
- **Categories:** smells
- **Problem:** `lay_out` and `part_repeat_places` each test the same three conditions; two came from a
  fix after the 2026-10-11 review.
- **Suggested fix:** one helper both call.

#### Mi8: [src/cli/inspect_psp.rs](../../../../src/cli/inspect_psp.rs#L269) — The empty rows' padding is hand-counted tab strings
- **Confidence:** High
- **Categories:** naming/smells, refactor_safety (nit)
- **Suggested fix:** derive the padding from `COLUMNS.len()`.

#### Mi9: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L1750) — "The review's finding M2/M3" names no review
- **Confidence:** High
- **Categories:** naming (clear-technical-writing)
- **Problem:** the same file already has an older "the review's M3" (line 4103) about a different
  finding.
- **Suggested fix:** name the review file at each citation.

#### Mi10: [src/locus_generation/mod.rs](../../../../src/locus_generation/mod.rs#L387) — `held_flank` can be set on a complete observation
- **Confidence:** Medium
- **Categories:** naming/idiomatic
- **Suggested fix:** move it into `ReadWitness::Partial`, or at least refuse it at the decoder (Mi1).

#### Mi11: [src/locus_generation/ssr.rs](../../../../src/locus_generation/ssr.rs#L2704) — `ObservationIdentity` is a five-element tuple
- **Confidence:** High
- **Categories:** idiomatic
- **Suggested fix:** a struct with named fields, deriving `Hash` and `Eq`.

#### Mi12: [src/calling/likelihood/ssr_emission.rs](../../../../src/calling/likelihood/ssr_emission.rs#L1724) — The ½/½ fallback for an unknown flank is silent
- **Confidence:** Medium
- **Categories:** errors
- **Problem:** every whole-hull partial read from a 2.0 psp is scored from both flanks, with no count
  in the run report. A cohort mixing 2.0 and 2.1 psps scores the same evidence two ways by sample.
- **Suggested fix:** count such reads in the run report.

#### Mi13: [doc/devel/ng/spec/str_bundle_loci.md](../../../../doc/devel/ng/spec/str_bundle_loci.md#L945) — "Both rows are over the same psps" is not true of the HG002 column
- **Confidence:** High
- **Categories:** extras (diff vs stated intent)
- **Problem:** HG002 is called from the BAM by each build; its step 4b cell came from an earlier run.
  The step 4c VCF also has one record more than step 4b (20,348 against 20,347), which the text does
  not say.
- **Suggested fix:** say so in §10.

### Nits

- `partial(..., border: HeldFlank)` keeps the deleted type's word; `held` reads better
  (`src/locus_generation/ssr.rs:1675`).
- `HeldFlank`'s two external spellings, the psp code and the TSV word, could be methods on the type.
- A leftover `eprintln!` in `a_read_stopping_before_a_run_scores_alleles_differing_there_alike`; two
  `#[allow(clippy::too_many_arguments)]` with no reason.
- The memo's `NaN` sentinel would recompute a computed `NaN` forever (none can arise today); and the
  held-flank tie-break test has two buckets, so a deleted `then_with` passes about half the time
  under hash-map seeding (float_portability).
- Stale "identity is `(bases, read_witness, read_group)`" comments: `witness.rs:85`, `mod.rs:292`,
  `ssr.rs:2721`, `:2735`, `parity.rs:833`.
- A comment quotes "4 parts in 10¹³" where the assertion checks below 10⁻⁹
  (`ssr_emission.rs:4542`).

## 7. Out of scope observations

- `HeldFlank` would sit better in `locus_generation/witness.rs`, beside `ReadWitness::from_left` and
  `from_right`. More broadly, the observation types shared by `psp`, `calling` and `run` live in a
  pipeline-stage module rather than a peer module (refactor_module reviewer).
- The flank-deletion decoder has the same leniency as Mi1.

## 8. Missing tests to add now

- `classify` — `a_partial_over_the_whole_hull_from_the_right_records_the_right_flank` (B1).
- `genotype_log_likelihood_row` — `the_row_hands_the_model_each_partials_held_flank` (M2).
- `bundle_censored_emission` — `a_partial_read_is_compared_at_the_rates_of_what_it_reached` (M1):
  bit-for-bit equal scores across counts of an unreached run under per-count rates, from both
  flanks; and the same read at the cell's mean rate moving with the far count.
- `bundle_censored_emission` — `a_two_run_partial_read_is_at_least_as_likely_as_the_complete_read` (M3).
- `gather_for_bundle` — `a_bundle_cell_carries_each_runs_rate_and_a_fixed_rate_no_candidate_moves` (M1).
- `decode_record_body` — duplicate index and a flank on a complete observation refused (Mi1).
- `check_this_build_writes` — refuses to extend a 2.0 file (Mi2).
- `inspect_psp` — `inspect_psp_writes_the_held_flank_column` (Mi3).

## 9. What's good

- The neutrality claim is pinned bit for bit, not to a tolerance, over alleles that differ beyond the
  read's reach from both flanks (`a_partial_read_scores_alike_against_alleles_that_differ_beyond_its_reach`).
- The psp bump is a genuine minor: verified from the HEAD code, a 2.0 reader takes minor 1, matches
  the 2.1 manifest as its own plus one unknown final field, skips it, and scores whole-hull partials
  ½/½ as it always did; the step 4b binary read all three 2.1 psp sets on the benchmark.
- The flank is recorded only where the witness cannot say it and only at bundles, so every segment
  observation and its key are unchanged — confirmed by 156,000 segment cells identical to main.
- Every production construction of `SequenceObservation` and `PartialObservation` is an exhaustive
  literal, so the compiler found every site for the new field.

### Quantitative claims in the diff

Every checkable figure was re-derived (errors/extras reviewer) and is correct: 184 / 1,030 / 50,776
whole-hull partials, about one in eight of the partial reads (1 in 7.7, 8.4, 8.5); their split
between flanks (90/94, 510/520, 25,562/25,214); both rows of the step 4c table; 59 of 1,898,
151 of 2,591, 452 of 6,005 and 1,336 of 20,348 records changed; the cost ratio 1.371, about 9.7 ms
more a bundle locus; 0.960567 equal to 13 significant digits across the five alleles; "one byte"
for the empty field. Not checkable from logs: "load average 55".

## 10. Commands to re-verify

- `scripts/dev.sh bash tmp/step4c_review/check.sh` — fmt on the touched files, clippy on every
  target, the touched modules' tests and the integration test, the release build.
- `scripts/dev.sh bash tmp/step4c_review/cost.sh 3 <bin> ...` — alternating one-thread CPU timings on
  the 500-bundle cost set, with a byte comparison of the VCFs.
- `scripts/dev.sh bash tmp/step4c_review/profile.sh <profiling-bin> <label>` — the cpu-clock profile.
- `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks scripts/dev.sh bash tmp/step4c_review/bench.sh`
  — the three benchmark runs and HG002 30×.
