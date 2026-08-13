# Fixes applied — ng parameter pre-pass, the STR path, A1+A2+A3

*2026-08-11, against [the review](ng_parameter_prepass_ssr_a1a2a3_2026-08-11.md) of steps A1+A2+A3
of [`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md). Eight agents, 46
mutations, 4 survivors, 0 Blockers, 4 Majors, 12 Minors.*

## Applied — 15

**M1 — `OffsetBucket` moved to its own file.** New `src/ng/parameter_estimation/ssr/offset_bucket.rs`
holds `OffsetBucket`, `bucket_of`, `OFFSET_HALF_RANGE` and `OFFSET_BUCKETS`; `mod.rs` re-exports
them. The three files that will index an entry's counts are now the bucket's **siblings** rather
than its descendants, so its private field constrains them: `OffsetBucket(200)` is `error[E0603]`
where it previously compiled and panicked at run time. The file's own module doc records why it is a
file.

**M2 — `SsrPeriod::try_new` takes a `usize`**, and `DomainError::SsrPeriod` carries one. The stored
field stays `u8`; the narrowing happens after the check, where it cannot lose anything. A test at
258 pins the value a `u8` parameter would have accepted as a dinucleotide.

**M3 — two compile-time assertions.** `ALLELE_OFFSET_LIMIT > OFFSET_HALF_RANGE` in `mod.rs`, where
both are in scope, and `OFFSET_HALF_RANGE > 0` beside the constant itself. Each fails as
`error[E0080]` rather than as a run-time surprise, which is the shape `generic/depth_bins.rs`
already uses.

**M4 — the two constants A3 specifies and this step had not written**: `MAX_LOCUS_READS = 12` and
`GUARD_SHARE_LIMIT = 0.10`, each with the warrant its architecture entry carries and a line naming
the milestone that consumes it (C and B respectively).

**Mi1 — the precedence trap.** `let lowest = -(reach_down.min(ALLELE_OFFSET_LIMIT));` with the
intermediate named `reach_down` and a comment stating what the near-identical misreading would
return: −6 at three reference repeats, a support reaching below an empty tract.

**Mi2, Mi3, Mi4, Mi5 — four tests the mutations asked for**: a proptest over `0..=1000` periods; the
allele support asserted as a whole vector at 0, 3, 5, 6 and 20 reference repeats; every one of the
256 `i8` offsets asserted in range and non-decreasing; and the offset's sign asserted in both
directions and at the origin.

**Mi6 — `a_stratum_cannot_be_built_at_a_period_of_zero` replaced** by
`a_stratum_carries_back_the_period_and_repeat_count_it_was_built_from`, which tests what this module
owns. The rejection it duplicated stays in `types.rs`, where the check is.

**Mi7 — `Display for Stratum` destructures**, so a third field is a compile error here rather than a
rendering that names two thirds of the key.

**Mi8–Mi11 — the four wrong numbers corrected** in the implementation report, with the correction
itself visible: the full-suite row now gives 3,386 across 12 binaries and says which line the
earlier draft had quoted; the clippy paragraph gives 8 errors, names all six lint kinds, and records
that the printed count is unstable because `-D warnings` aborts at the first failing example.

**Mi12 — the module doc's definition of a stratum**, which said three parts where the type has two.
It now says a stratum is the period and the repeat count, and that the read group is the other half
of the key the fits are filed under.

**Nit — `#[must_use]`** on the two free functions and the accessors, matching the sibling path.

**Out-of-scope item taken anyway:** the new doc comments no longer cite `research note §6.8` and
`§6.4.1`, which do not resolve (Open question 1). They cite the spec sections that actually hold the
figures. **The design documents still carry the unresolvable citations** — correcting those is the
owner's call at Checkpoint A, because the question is not the citation but whether the numbers have
a source at all.

## Not applied — 4, each with its reason

- **`RepeatCount` → `ReferenceRepeatCount`** (naming, Minor). The architecture names the type
  `RepeatCount` and its doc comment already says "the reference tract's count, never the sample's".
  Deferred to Milestone C, which is where a read's own count first appears beside it and where the
  ambiguity the finding predicts would actually bite.
- **`allele_support` returning an iterator instead of a `Vec`** (idiomatic, Minor). It is called
  once per stratum per fit — a few hundred times a sample — against a search that re-walks a
  stratum's whole table per candidate. Milestone D indexes the support pairwise for genotypes, which
  a `Vec` serves directly. Recorded rather than taken.
- **`bucket_of` → `OffsetBucket::recording`** (naming, Minor). The plan names the function
  `bucket_of` in its A3 step and in the test it specifies; keeping the plan's name keeps the step
  checkable against the plan.
- **`#[cfg(test)]` blocks in the three scaffold files** (reliability, Nit). An empty test module
  asserts nothing. Already recorded as a deviation in the implementation report.

## Validation after the fixes

All in the container:

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --all-features` | **3,391 passed, 0 failed, 10 ignored** across 12 binaries (3,311 the library's) |
| `cargo test --lib parameter_estimation::ssr` | 14 passed |
| `cargo test --lib ng::types` | 34 passed |

**19 tests new in this step** — 8 in `ssr/mod.rs`, 6 in `ssr/offset_bucket.rs`, 5 in `types.rs`
(four named, one proptest). Five constants: three in `mod.rs`, two in `offset_bucket.rs`. Every
number in this paragraph was counted from the tree rather than recalled — `grep -c '#\[test\]'` and
`grep -c 'pub const'` on the files, and the test totals from the run above.
