# Code review — ng cohort merge, B3: `CohortObservation` and per-sample support

*2026-08-17. Six category checklists, each in its own isolated git worktree detached at the
commit under review. Per-category audit trail in
`tmp/review_2026-08-17_ng-cohort-merge-b3/`.*

## 1. Scope

- **What:** the working-tree diff of step B3, as commit `67b55cb2` (parent: the B2 commit,
  branch `ng-cohort-merge`). One file changed, `src/ng/run/cohort_merge/build.rs`, +714/−10.
- **In scope:** `CohortObservation`, `SampleSupport`, `AlleleSupport`, `AlleleTable::assemble`,
  the division (`share_of_one_read`, `AlleleSupportTally`, `round_to_u32`/`round_to_u64`), the
  derivation's new callback argument, and the tests.
- **Out of scope:** `close.rs`, `mod.rs`; `src/var_calling/`, `src/pileup/`, `freebayes/`
  (read-only references); the organiser and the parallel arrangement.
- **Categories dispatched:** reliability (the mutation pass), errors (floating point enters
  the module here), idiomatic + naming (this is the shape every later step consumes),
  refactor_safety (`AlleleTable::over` stopped being the builder), smells, extras (intent,
  hot path, stable output).

## 2. Verdict

**Approve with changes** — 2 Blockers, 4 Majors, 14 Minors, 12 Nits. Both Blockers were
tests that could not fail. One Major is a measurement that needs the owner (§5).

## 3. Execution status

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean. (`--all-targets` red on this
  branch, identically before this work.)
- `cargo test --lib ng::run::cohort_merge` — 94 passed at review time.
- `cargo test --lib` — 3,717 passed, 0 failed, 11 ignored.

## 4. The sign question, checked three ways

The step's most consequential claim is that production's `min` over the constituents' mean
`q_sum` picks the **best** piece of a compound and contradicts production's own plan, so ng
takes the `max`. Independently confirmed:

- `q_sum` is `Σ ln P(error)` (`locus_generation/mod.rs`, `read/prepared_read.rs` — `mq_log_err`
  is "`ln(P_err)`"), so it is ≤ 0 and a weaker read is a **larger** number: 1 error in 400 is
  −6.0, 1 in 3 is −1.0.
- **ng's own pileup already reduces two error terms with `max` for the same reason** —
  `ln_bq_for_read(...).max(mq_log_err)` (`pileup/open_record.rs`), the read's effective error
  being the worse of its base-quality and mapping terms. B3 applies that reduction across
  records instead of across terms.
- "Quality cannot exceed any constituent's" translates to `q_sum ≥ every constituent's`,
  which is `max`.

## 5. The Majors

| # | Finding | Resolution |
|---|---|---|
| M1 | **The dense support row is the module's memory at the top of the cohort range.** Measured: a locus where every sample shows a distinct allele costs **614 MB for one observation at 4,000 samples**, against 41 MB at 1,000 — the `samples × alleles` matrix. B3 is otherwise cheap (1.6 → 3.1 µs at 63 samples × 3 reads; +11% on the cross-record shape). Spec §8 prices a survivor with no sample-count factor at all. | **Raised at Checkpoint B, not decided.** Sparse row, a cap on alleles, or pricing it as it stands is a design choice. Recorded at the field with the measurement. |
| M2 | **A producer-guarantee check without the psp obligation, and the same defect loud on one branch and silent on the other**: `share_of_one_read` asserted a sequence's read count is non-zero, while the one-record branch happily interned a zero-read sequence as an allele carrying a quality nobody measured. | **Applied.** Both branches now skip a sequence no read is behind; the assertion is a documented backstop and carries the `RunError` obligation its three siblings carry. |
| M3 | **The read-group axis does not survive into a cohort observation** — `SequenceObservation` keeps it because "a per-chemistry model needs the allele × group cross with its quality moments", and one row per allele pools it. The reviewer confirms arch §4's own sketch pooled it too, so B3 made the loss visible rather than causing it. | **Recorded at the type**, with what wants it (the pre-pass fits ε per read group, from the census rather than from these observations) and the shape that would restore it. For the owner. |
| M4 | **`placed_left` is counted against each *record's* position, not the locus's**, so pooling it across records mixes two questions — while the doc claimed the locus's anchor and justified the pooling with "properties of the read, identical at each sighting", which is true of MAPQ and strand and false of placement. | **Applied**: the field says what it is, and the division's justification is split into the two cases with the approximation's bound named (the locus width). |

## 6. The two Blockers — both tests that could not fail

**B1 — the quality rule's fixtures could not tell `max` from "keep the last sighting".** All
three put the weakest sighting last, so `Some(_) => mean_quality` left 94 tests green.

**B2 — which sequence a sighting names was never exercised.** Replacing
`observations[sighting.sequence]` with `observations[0]` left 94 tests green, because no
fixture had a record with two sequences on the composed path — **so the case the whole
division exists for, one observation's reads splitting onto two alleles, had no test at all**.

Both are closed by one fixture,
`one_observations_reads_split_across_two_alleles_and_each_takes_its_own_share`: a record with
two sequences, three reads taking two paths, the weakest sighting **first** for two of them
and **last** for the third. Both mutations were re-run against it and now fail.

## 7. Other findings applied

- **Two silent fallbacks for states the code's own assertion rules out** — `unwrap_or(0.0)`
  on the quality (and `0.0` is not neutral: `ln P(error) = 0` is the worst quality
  expressible) and a guarded division. Both removed.
- **The five-field support list was restated in eight places**, only one of which fails to
  compile if a sixth is added. Now one type (`SupportSums`) written once for gathering and
  once for rounding, with the adds destructured so a new field is a compile error.
- **A dead `resize`**, proved redundant with the level-up ten lines below.
- **`round_to_u32`'s doc was wrong about the language**: `as` has saturated since Rust 1.45,
  so these helpers are not a repair. Doc corrected and the boundaries pinned by a test.
- **`num_reads` was labelled "Exact." while built with `saturating_add`**;
  `reads_composed_across_records` omitted "Saturating."; `reads_removed_as_evidence` omitted
  the mate-overlap removal case. All corrected.
- **`ShownBy` was a preposition** (the project's types-are-nouns rule) with variants at
  different grains → `AlleleBacking::{OneSequence, OneRead}`, and **the records now travel
  inside the `OneRead` variant**, since a sighting is a pair of indices that means nothing
  without the slice it indexes — the same defect B1's review found in the projection.
- **The per-allele tally is hoisted and refilled** rather than allocated per sample.
- **`CohortObservation::over` destructures the table** rather than reading one field, so
  anything the table gains has to be answered for.
- **A wrong mechanism in my own comment**: the zero-read guard claimed the division "would
  come back as an infinity and poison every score, silently". Measured with the guard
  removed, the `-inf` is **discarded** by the `max` and the allele reports a plausible
  quality — worse than the story the comment told, and now the story it tells.

## 8. Declined, with reasons

- **Deleting `round_to_u32`/`round_to_u64` in favour of bare `as`.** They now say what they
  do, have a test that pins `NaN`, both saturation ends and the half-way rounding, and cost
  nothing.
- **`AlleleSupport`/`SampleSupport` renamed to fix the grammatical clash** (an allele is
  supported; a sample supports). `SampleSupport` is the name arch §4 fixes — the owner's to
  settle, with the arch.
- **A private field on `CohortObservation` to protect the parallel-vector invariant.** It
  changes the shape arch §4 declares.
- **Sparse or capped support rows** — M1, the owner's.

## 9. What's good

- The refactor was proved rather than argued, twice: the B2 `AlleleTable::over` and
  `alleles_of_sample` restored under second names and compared against today's over **844
  generated loci** (14 sample shapes, one to three samples, two coordinate bases including
  one whose locus ends exactly at `u64::MAX`) — alleles, order, emission sequence and removal
  count identical, with the differentials themselves mutation-tested so that passing means
  something.
- **Determinism proved three ways** rather than asserted: rebuilding with every sequence's
  ids and every record's sequence list reversed gives a byte-identical dump at 17 digits; no
  `f64` crosses samples; and six separate processes agree on the allele table's hash, which
  matters because `ahash` seeds per process.
- The mutation pass ran 21 mutations and reported all three numbers — 21 run, 6 survived, 1
  changed no behaviour — and paired every survivor with a fixture proving the mutant answers
  differently.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh cargo test --lib
```
</content>
