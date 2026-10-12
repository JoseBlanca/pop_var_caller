# Code Review: str_partial_read_phase
**Date:** 2026-10-12
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** uncommitted diff on branch `str-partial-phase` (from main a87ce7b6): a partial read at an STR segment is compared with the end of the candidate it held
**Status:** Request-changes (one Blocker, two Majors; all three fixed in the fix run that followed — see `fixes_applied_str_partial_read_phase_2026-10-12.md`)

---

## 1. Scope

- **What was reviewed:** the working-tree diff of `/Users/jose/devel/pop_var_caller-str-partial-phase`
  against a87ce7b6, and the untracked `scripts/ng_str_accuracy/hg002_two_binaries.sh`.
- **Reviewed against:** a87ce7b6 plus the uncommitted change, as it stood before the fix run.
- **In-scope files:** `src/calling/likelihood/ssr_emission.rs`, `src/calling/likelihood/ssr.rs`,
  `src/calling/likelihood/mod.rs`, `src/locus_generation/mod.rs`, `src/locus_generation/witness.rs`,
  `src/calling/evidence_shaping.rs`, `src/calling/mod.rs`, the three `src/calling/inference/` files
  (test call sites only), `src/cli/cross_platform_digests.rs`, the dated notes in
  `doc/devel/ng/spec/read_likelihoods.md` §5.2 and `locus_generation_ssr.md` §3, and the
  implementation report `doc/devel/reports/implementations/str_partial_read_phase_2026-10-12.md`.
- **Deliberately out of scope:** the STR generator's classification of partial reads (unchanged);
  the complete-read emission (unchanged, pinned by a new test).
- **Categories dispatched.** **Owner-required deviation from the skill: no git worktrees.** Every
  category reviewer read the code in this worktree and neither built nor mutated it; the probes that
  decide the Blocker were run by the orchestrator afterwards.
  - `reliability`, `errors`, `naming`, `idiomatic`, `refactor_safety`, `smells` — always.
  - `module_structure` — the change spans two top-level modules (`locus_generation`, `calling`); it
    also carried the comparison with branch `str-bundles`.
  - `float_portability` — the change adds float arithmetic that reaches calls and re-records two
    cross-platform checksums.
  - `extras` — a PR-shaped change: does the diff match its report and spec notes.
  - `defaults` — a new literal weight (one half) acts as a default.
  - Not dispatched: `unsafe_concurrency` (no `unsafe`, no threads in scope), `tooling` (no
    `Cargo.toml` change).

## 2. Verdict

Request-changes.

## 3. Execution status

Run by the orchestrator in the dev container, after the read-only fan-out:

- `cargo test --release --no-fail-fast --lib -- a_right_held_read_of_a_candidate_ending_mid_copy`
  on the reviewed code → exit 101: `8 bases: right 0.00337285823656753, mirror from the left
  0.003541442065278368`. This is the Blocker, reproduced.
- `cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run::callers
  cli::cross_platform_digests` (after the fixes) → exit 0, `1674 passed; 0 failed; 3 ignored`.
- `cargo fmt --check` on the two touched source files → exit 0 after formatting.
- `cargo clippy --release --all-targets -- -D warnings` → exit 101, but on nothing this change
  touches: `type_complexity` at `src/locus_generation/ssr.rs:2117` and `single_match` at
  `src/run/edge_deletion_gathering.rs:201`, both pre-existing; one `type_complexity` in a new test
  was fixed.

Not run: `cargo doc`, `cargo audit` (no dependency or public-doc change of note).
Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. **The one-half weight for a whole-tract witness** (affects Mi4, Mi8). A read whose witness
   covers the whole tract is scored as an even mixture of "held the left flank" and "held the right
   flank". As a generative model this is right when nothing else is known, and on a pure tract it
   costs only a factor of one half that every candidate shares. It disappears once the generator
   records the side (str-bundles' `held_flank`, psp 2.1). Accepted as is; recorded here so the
   owner sees it was looked at.
2. **Partial reads against alleles that delete flank bases** (`str_edge_deletion_alleles.md`). A
   partial read carries no flank deletion, and the censored term compares it with the allele's
   spelling as stored. Reading `respelled_over`, an allele's bases are its carrier's sequence
   between the flank bases it keeps, so comparing a right-held read with the spelling's end looks
   correct; not measured, and the implementation report says so (§7).
3. **The greedy segmentation, mirrored, splits interrupted tracts from the right.** For a
   candidate whose runs split differently from either end, a right-held read is scored over
   different placements than a left-held or complete read of the same candidate. The spec note
   states this ("as found from the right"); it is a convention, not a defect.

## 5. Top 3 priorities

1. **B1** — a right-held read of a candidate that ends part-way through a copy is mirrored out of
   phase; the change's own claim (an exact mirror image) fails there, and no test reached it.
2. **M2** — a partial witness touching neither border of the tract is an impossible state, and it
   was scored half from each side without a word; it is also the symptom of a tract length on the
   wrong axis.
3. **M1** — the one interrupted-tract mirror test is a palindrome at period 1, so it cannot see a
   missing reversal of the candidate or the motif.

## 6. Findings

### Blocker

- **B1: `src/calling/likelihood/ssr_emission.rs:465-472` — A candidate that ends part-way through a copy is mirrored out of phase**
- **Categories:** reliability (as a test gap), orchestrator (the mechanism, reproduced)
- **Confidence:** High
- **Assumptions:** none.
- **Problem:** The right-held branch writes the motif plainly backwards. A tract written backwards
  starts where it used to end, so its motif is the copy it **ended** on, written backwards; that is
  the motif reversed only when the candidate ends on a whole copy. `CACACAC` under `CA` ends on
  `AC`; its mirror image `CACACAC` is in phase with `CA`, but was handed `AC`. `resize_at_the_end`
  tiles by offset from the start, so a lengthening of part of a repeat appends `C` where the mirror
  image of the left-held rule appends `A`. Measured on that candidate, 8-base right-held read: 3.373
  in 1,000 against the mirror image's 3.541 in 1,000. Reference tracts are trimmed to whole copies
  (`repeat_catalog.md`, the "trimmed start, end" row), so only candidates proposed by the reads are
  affected, and only on their part-repeat terms.
- **Why it matters:** wrong likelihoods on the path the change exists to fix, and the doc comment's
  claim — "every length change lands as the mirror image" — false for this class. Every new test
  built its candidates with `a_tract`, whole copies only, so none could see it.
- **Suggested fix:** derive the mirrored motif from the candidate's last `period` bases when they
  are the motif in some phase, else the motif reversed; and test `score(Right, r, c, m) ==
  score(Left, reverse(r), reverse(c), m')` with `m'` written out by hand, on `CACACAC` and on
  interrupted candidates that are not palindromes.

### Major

- **M1: `src/calling/likelihood/ssr_emission.rs` test `mirror_image_reads_of_a_symmetric_interrupted_candidate_score_alike` — cannot fail for what it claims**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** its doc says it pins whole-repeat slips placed run by run and a part-repeat change at
  the far end. The fixture is period 1, so no part-repeat change exists; `AAAATAAAA` is a
  palindrome, so a missing reversal of the candidate gives the same bytes; and `A` is its own
  reverse. An implementation that reversed only the read passes.
- **Why it matters:** it reads as coverage of the interrupted route under the mirror, which is where
  B1 lives.
- **Suggested fix:** an interrupted fixture at period ≥ 2 that is not a palindrome, with the
  mirror-image motif written out by hand (folded into B1's test).

- **M2: `src/locus_generation/mod.rs:399-412` and `src/calling/likelihood/ssr.rs:861-876` — a witness touching neither border is silently scored half from each side; the tract length is unchecked**
- **Categories:** errors, refactor_safety, reliability, extras (convergent)
- **Confidence:** High (the state is impossible from the generator); Medium (that the two lengths
  can diverge today).
- **Problem:** `HeldFlank::of_witness` answers `None` for a whole-tract witness, a complete read
  and a witness touching neither border alike, and the row scores every `None` as an even mixture.
  The STR generator never mints a neither-border witness (`from_left`/`from_right` always touch one
  border, `src/locus_generation/ssr.rs:1461-1463`). It appears only when the tract length handed to
  calling (`LocusLen::of_region(region)`, `evidence_shaping.rs:464`, the merged record's region) is
  longer than the one the witness was laid on (`segment.tract_len()`, `ssr.rs:1447`). They agree
  today: the record region is the segment, and `the_mints_locus_length_is_the_emitted_regions`
  pins it for one sample. But nothing at calling checks it, and a divergence would silently turn
  every right-held read back into a blend — undoing this change without a failing test. The report,
  the spec note and the code comment all describe the `None` branch as the whole-tract case only.
- **Why it matters:** an error that cannot pass silently is passing silently, on the path that
  decides which end of the candidate a read is compared with.
- **Suggested fix:** in the row, refuse a partial whose side is `None` unless its witness is flush
  with both borders, naming the cause; a `should_panic` test laying a right-held witness on 8
  positions and scoring it against 10. (A shorter calling length — a left-held read becoming
  whole-tract — is not caught by this; see Residual risk in the fix report.)

### Minor

- **Mi1: `src/locus_generation/mod.rs:399` — `HeldFlank::of_witness` has no direct test.** Complete,
  neither-border and holed both-border witnesses are reached by no test. Add
  `of_witness_answers_none_for_complete_interior_and_holed_witnesses`. (reliability, extras)
- **Mi2: `src/locus_generation/mod.rs:386` — `HeldFlank` duplicates the producer's private
  `AnchoredBorder` (`src/locus_generation/ssr.rs:1384`).** The two halves of one round trip under
  two names; replace `AnchoredBorder` with `HeldFlank` and add a round-trip test. Best done after
  the str-bundles merge, which touches the same enum. (smells, module_structure)
- **Mi3: naming.** The closure `from` in `fill_emissions` (`ssr.rs:864`) reads as `From::from`;
  `witnessed` is a bare participle and the rename stopped short of the counting mock and
  `score_censored_from`; `tract_len` names the value `SampleLocusObservations::locus_len()` calls
  `locus_len`; test bindings `as_main_read_it` and `main` name the branch and invert their meaning
  once merged. (naming)
- **Mi4: `src/calling/likelihood/ssr.rs:875` — the weight one half is an unnamed literal.**
  (defaults)
- **Mi5: `src/calling/likelihood/ssr.rs:757-880` — `fill_emissions` passes 120 lines;** extract the
  partial-read loop. (smells)
- **Mi6: `src/calling/likelihood/mod.rs:1575` — `TEST_TRACT_LEN` of 1,000.** It does not hide any
  existing right-held coverage (no fixture passing through `SsrSampleEvidence` built one — the only
  `from_right` under `src/calling/` is the new test, which uses its own length), but a future
  right-flush fixture would silently read as neither border. After M2's fix that fixture panics
  instead, which is the better failure. (refactor_safety, idiomatic, reliability)
- **Mi7: `src/cli/cross_platform_digests.rs` notes for `WITH_ERRORS_CALLS_MD5`,
  `OF_21_SAMPLES_CALLS_MD5` — the isolating check is in the implementation report (§6: scoring
  every partial from the left again restores both old checksums) but not in the notes,** as earlier
  notes in the file have it. (float_portability, extras)
- **Mi8: implementation report §4 — the explanation of why homopolymers moved is half wrong.** It
  says homopolymers move because an interrupted run of one letter, *or a whole-tract read now scored
  from both sides*, is not its own mirror image. For a pure run of one letter the two sides are
  equal: a reversed read has the same mismatch count against every candidate. Traced on the 50
  samples, 5× run: 20 period-1 cells changed verdict and 18 of them have a truth of
  `substitution_inside` (an interrupted run); at 30× one cell moved in each run, both
  `substitution_inside`. So period-1 cells should move, and they move for the first reason only.
  (orchestrator)
- **Mi9: `scripts/ng_str_accuracy/hg002_two_binaries.sh` — `|| true` hides a failure of the
  verdict-flip scorer, the dump binary is looked for only under `target-container/` (CLAUDE.md asks
  scripts to check `target/` too and take the newer), and its environment defaults are neither
  printed nor part of the cache key.** The script is byte-identical to the copy str-bundles
  committed in 199541ac. (errors, defaults, reliability)
- **Mi10: `src/calling/likelihood/ssr.rs:875` — no test scores both sides through one reused
  scratch.** Every buffer is cleared before use today; a later missing `clear()` would not fail a
  test. (reliability)

### Nits

- The bit-equality assertions in the mirror tests hold because the fixtures' summed terms are
  identical sequences; a future mirror test with asymmetric errors should compare within a few
  units in the last place (float_portability).
- Change history in the trait's doc comment belongs in the spec box it cites (smells).
- `HeldFlank` derives `Hash`, `PartialOrd`, `Ord` that nothing uses; import placement of
  `TEST_TRACT_LEN` above `std` imports; the `of_witness` match could spell its two `None` rows
  (idiomatic).
- `censored_from_the_tract_start` holds two motifs in scope, one possibly reversed, kept apart by
  a comment (smells).

## 7. Out of scope observations

- `src/calling/likelihood/ssr_emission.rs`, `resize_at_the_end` — its doc says lengthening
  "continues the motif's tiling from the phase the tract ended on", but it tiles by offset from the
  tract's start. The two differ after an interruption whose length is not a whole number of repeats
  (`CACATCACA` ends on `CA`; offset 9 tiles `A` next, where the tract's own phase gives `C`). This
  affects complete and left-held reads too, predates this change, and the mirror inherits it
  symmetrically. Follow-up: a separate change, measured.
- `src/locus_generation/ssr.rs:2117` and `src/run/edge_deletion_gathering.rs:201` fail
  `clippy -D warnings` on main's code.

### Conflict with branch `str-bundles`, and the merge order

Both branches add the same `HeldFlank` enum, word for word, at the same place in
`src/locus_generation/mod.rs`; str-bundles also adds `SequenceObservation::held_flank:
Option<HeldFlank>` (psp 2.1) in that hunk. The real conflicts are elsewhere:

- **The trait.** Here `censored_emission` takes a definite `HeldFlank` and the row decides; on
  str-bundles it takes the witness positions and the recorded `Option<HeldFlank>` and the model
  decides. Only one can stand.
- **Where the even mixture lives.** Here in the row; there inside the bundle emission. A naive
  merge blends a bundle read twice or not at all.
- **Two functions decide the side.** `HeldFlank::of_witness` here; a private `held_flank_of` there,
  which falls back to the recorded side.
- **Two length sources.** The record region here; `bundle.reference_hull_len` there.
- **Stale documentation.** str-bundles says a segment reads every partial from the tract's start;
  false once this lands.
- **The segment's censored body** moved into `censored_from_the_tract_start` here; kept inline
  there. About 25 small test-site conflicts besides.

**Recommendation: merge `str-partial-phase` to main first, then rebase str-bundles onto it.** This
branch is small, measured and independent; str-bundles is large (136 files), still has
uncommitted work, and its psp 2.1 field is what removes this branch's one assumption. On the
rebase: keep one `HeldFlank` and one side function `of_witness(witness, locus_len, recorded)`;
one trait signature taking a definite side, with the even mixture in the row; one source for the
locus length; and let the generator record `held_flank` for segments too, so a whole-tract witness
at a segment is no longer blended. Add a right-held arm to str-bundles'
`a_one_run_bundles_partial_read_scores_as_the_segments` to show the two mirror mechanisms agree.

## 8. Missing tests to add now

- `ssr_emission`: `a_right_held_read_of_a_candidate_ending_mid_copy_mirrors_the_left_held_score`
  — candidates `CACACAC`, `CACACAGCACA`, `CAGCAGCATCAGCA`; would catch B1 and M1. (Added.)
- `ssr`: `a_partial_witness_touching_neither_border_is_refused` — would catch M2. (Added.)
- `locus_generation`: `of_witness_answers_none_for_complete_interior_and_holed_witnesses` (Mi1).
- `ssr_emission`: `censored_emission_is_the_same_with_a_scratch_reused_across_sides` (Mi10).
- `evidence_shaping`: `shape_ssr_locus_gives_every_sample_the_regions_length` — a right-held
  partial minted on the region's own length comes back `HeldFlank::Right` (M2's other half).

## 9. What's good

- Turning a right-held read into a left-held one by reversal keeps one scoring routine for both
  sides, so the stutter model, the placements and the outlier floor cannot drift between them
  (`ssr_emission.rs`, `censored_emission`).
- The scratch buffers are taken and put back (`mem::take` of `mirrored`), so the mirror allocates
  nothing after warm-up.
- The digest re-record was isolated by restoring the old rule and getting the old checksums back
  (implementation report §6).
- The row test with a stub model answering one value per side (`BySide`) pins the dispatch without
  depending on the emission's numbers.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run::callers cli::cross_platform_digests`
- `scripts/dev.sh cargo test --release --lib -- ending_mid_copy touching_neither_border`
- `scripts/dev.sh cargo clippy --release --all-targets -- -D warnings` (expect only the two
  pre-existing errors above)
