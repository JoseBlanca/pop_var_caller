# Code Review: ng parameter pre-pass, the STR path — A1+A2+A3
**Date:** 2026-08-11
**Reviewer:** rust-code-review skill (orchestrator), 8 sub-agents in isolated worktrees
**Scope:** the working-tree diff of steps A1+A2+A3 of `doc/devel/ng/impl_plan/parameter_prepass_ssr.md`
**Status:** Request-changes → resolved (see the fixes-applied report)

---

## 1. Scope

- **What:** the uncommitted working-tree diff of one implementation step — the `ssr/` module tree,
  `SsrPeriod` in the shared vocabulary, the stratum, the offset vocabulary and the two widths.
- **Reviewed against:** `b96fdcd0` plus `tmp/review_2026-08-11_ng-prepass-ssr-a1a2a3/step.patch`.
  Every agent verified it was reviewing that tree before starting; all eight passed the check.
- **In scope:** `src/ng/parameter_estimation/ssr/{mod,locus_offsets,stratum_table,slippage}.rs`,
  `src/ng/parameter_estimation/mod.rs`, `src/ng/types.rs`.
- **Out of scope:** the `doc/devel/ng/**.md` and `PROJECT_STATUS.md` prose in the same patch;
  `examples/`, whose clippy failures are pre-existing and were verified so.
- **Categories dispatched (8):** reliability, errors, naming, defaults, module_structure,
  idiomatic + smells (one agent, two sections), refactor_safety, and the skill's step 8a
  quantitative-claim check as its own agent.

## 2. Verdict

**Request-changes**, on one finding with a run-time consequence and four convergent ones about
constants and their relations. No Blockers. All resolved before commit.

## 3. Execution status

Run by the orchestrator in the container, quoted to the agents rather than re-run by each:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 101 | **8 errors, all in `examples/`, all pre-existing** — reproduced at `b96fdcd0` with the patch **not** applied |
| `cargo test --lib --bins --tests --all-features` | 0 | 3,386 passed, 0 failed, 10 ignored, across 12 binaries |
| `cargo doc --no-deps --lib` | 101 | 13 unresolved intra-doc links, all pre-existing |

Findings labelled "Needs verification": **0**. Every finding below was demonstrated by a mutation
that compiled, or by a command whose output the agent quoted.

**Mutation totals across the agents that ran them:** 28 + 10 + 8 = **46 mutants, 4 survived, 0
changed-no-behaviour**. Each survivor was proven to change behaviour by a probe test that failed
under the mutant and passed on unmutated code.

## 4. Open questions and assumptions

1. **Do the STR design documents' load-bearing table-size and offset-distribution figures have a
   source?** (affects nothing in this diff's code, everything in what it cites — see §7.) The
   spec, the architecture and the plan attribute "0.43 entries a locus", "12,727 entries for 29,811
   loci", "70,305 entries over 1.73 M tomato loci", "88.9% at the reference" and "the end buckets
   take 0.89% of reads" to `research/parameter_estimator_experiments_2026-08-06.md` §6.8 — **which
   is titled "What is not measured" and says the opposite**: *"Depth above about 12 reads a locus …
   the worlds here are exact to 12 reads. HG002 at 300× needs a coarsening of the per-locus cell
   that nothing here has priced."* There is no §6.4.1 either, which four citations name. The figures
   appear nowhere else in `doc/`. **For the owner at Checkpoint A.**
2. **Should `RepeatCount` be `ReferenceRepeatCount`?** The naming agent argues the type is silent
   about whose repeats it counts, when "the reference tract's, never the sample's" is exactly what
   makes strata comparable across samples. Kept as the architecture names it; recorded here because
   Milestone C is where a read's own count first appears beside it.

## 5. Top 3 priorities

1. **B/M1 — `OffsetBucket`'s private field does not constrain the three modules that will use it.**
   Rust privacy reaches descendants, so `OffsetBucket(200)` compiled inside `ssr/` and panicked with
   an out-of-bounds index at run time. (refactor_safety, convergent with naming and module_structure)
2. **M2 — `SsrPeriod::try_new` took a `u8` while every producer of a period holds a `usize`**, so
   the first caller writes `as u8` and 258 validates as a dinucleotide. (errors)
3. **M3 — the two widths are interchangeable `i8`s and their "wider than" relation was prose only.**
   Raising `OFFSET_HALF_RANGE` past `ALLELE_OFFSET_LIMIT` compiled clean and inverts the property
   that lets an end bucket be explained by an allele rather than a far slip. (defaults, convergent
   with refactor_safety, idiomatic and reliability)

## 6. Findings

### Major

**M1: `src/ng/parameter_estimation/ssr/mod.rs` — an `OffsetBucket` outside its documented range is
constructible by exactly the modules that will index with it.**
**Categories:** refactor_safety (demonstrated), naming, module_structure. **Confidence:** High.
The type's field was made private *because* an out-of-range bucket would index past an entry's
counts. Rust privacy is module-and-descendants, and `locus_offsets`, `stratum_table` and `slippage`
are declared as children of `ssr` — so all three sit inside the defining scope. Writing
`OffsetBucket(200)` in `stratum_table.rs` compiled and panicked at run time with
`index out of bounds: the len is 9 but the index is 200`; the identical line one module up is
`error[E0603]`. **Fix:** move `OffsetBucket` and `bucket_of` into `ssr/offset_bucket.rs` and
re-export, making the four workhorse files siblings rather than descendants.

**M2: `src/ng/types.rs` — `SsrPeriod::try_new(bases: u8)` is narrower than every caller that will
exist.** **Category:** errors. **Confidence:** High.
No such caller exists yet; A4, A5, B and C add them, and every producer of a period in the crate
holds a `usize` (`Motif::period`, a tract length divided by one). A `u8` parameter makes
`try_new(n as u8)` the natural call, under which `258 as u8 == 2` passes validation as a
dinucleotide. **Fix:** widen the parameter and the `DomainError` payload to `usize`; the stored
field stays `u8`.

**M3: `src/ng/parameter_estimation/ssr/mod.rs` — the relation between the two widths is stated in
prose and enforced by nothing.** **Categories:** defaults (demonstrated), refactor_safety,
idiomatic, reliability. **Confidence:** High.
`OFFSET_HALF_RANGE` and `ALLELE_OFFSET_LIMIT` are both bare `i8` in the same units, so each compiles
where the other belongs — both swaps were compiled. More seriously, setting `OFFSET_HALF_RANGE = 6`
compiles clean and inverts the property the design turns on; the spec prices the failure it opens at
+499% of the slippage level in the worst measured row. A negative half-range also compiles, giving
`OFFSET_BUCKETS ≈ 1.8 × 10¹⁹`. **Fix:** two `const _: () = assert!(…)`, the shape
`generic/depth_bins.rs` already uses for this class.

**M4: the step is two constants short of what A3 specifies.** **Categories:** module_structure,
reliability (both as cross-category notes). **Confidence:** High.
`doc/devel/ng/impl_plan/parameter_prepass_ssr.md` A3 lists `MAX_LOCUS_READS = 12` and
`GUARD_SHARE_LIMIT = 0.10`; neither appeared anywhere under `src/`. **Fix:** add both with the
doc comments their architecture entries carry.

### Minor

**Mi1: `allele_support`'s clip is correct only by operator precedence.** Flagged independently by
**five** agents — reliability, errors, idiomatic/smells, refactor_safety and the numbers check.
`-ALLELE_OFFSET_LIMIT.min(reach_down)` parses as `-(min(…))`, which is right, but reads as
`(-6).min(reach_down)`, which at 3 reference repeats gives −6: a support reaching three copies below
an empty tract. Two tests kill that mutant today; the expression still invites the edit. **Fix:**
name the intermediate and parenthesise.

**Mi2: `SsrPeriod`'s rejection was checked at 3 of the values it must reject.** **Category:**
reliability, demonstrated. Rewriting the bound as `(usize::from(bases) & 0x07) > MAX_MOTIF_LEN`
admits period 8 while 0, 7 and 255 all still reject — mutant survived the whole suite. Ten lines
above it in the same file, `Ploidy`'s proptest doc says exactly this: *"Named `every` and checking
three of 255 was the gap."* **Fix:** a proptest over the domain.

**Mi3: `allele_support`'s interior was unpinned.** **Category:** reliability, demonstrated. The test
asserted length, first and last, so a support of `[-3, -3, -1, 0, …, 6]` survived. A set with a
duplicate and a hole has the right cardinality, so the genotype count still comes out right while
the fit never considers one allele length and double-counts another. **Fix:** assert the whole
vector, and add `RepeatCount(5)` — the clip's last bite.

**Mi4: `bucket_of`'s totality was documented over the whole `i8` domain and tested at 17 points.**
**Category:** reliability. No mutant escaped through this gap, so it is filed on the invariant
rather than on demonstrated escape — but the index subscripts a fixed-size array from Milestone B
on, and 256 iterations settle it permanently. **Fix:** the exhaustive loop, also asserting
monotonicity.

**Mi5: `WholeRepeatOffset`'s "signed, always" had no test.** **Category:** reliability,
demonstrated — `{:+}` → `{}` left the suite green, while the module's other three `Display` impls
are all covered.

**Mi6: `a_stratum_cannot_be_built_at_a_period_of_zero` asserted nothing about a stratum**, and
duplicated a `types.rs` test line for line. **Category:** reliability. No change to `Stratum` can
make it fail, so the property its name claims is enforced by the type system and observed by
nothing. **Fix:** test what the module owns — that a stratum carries back the pair it was built
from.

**Mi7: `Display for Stratum` read its fields rather than destructuring**, so a third field would
leave the rendering silently two-thirds complete. **Category:** refactor_safety, demonstrated
(exactly one compile error, and not in the impl).

**Mi8–Mi11: four wrong numbers, every one the author's own claim about the author's own gate run.**
**Category:** step 8a. The full-suite row quoted one binary's line (3,306) as the whole command's
(3,386 across 12); the pre-existing clippy failure is 8 errors, not 7, and the printed count is
unstable because `-D warnings` aborts at the first failing example; 2 of the 8 lints are the kinds
named as "all of them"; and "the four constants" introduced a list of three. **All 22 figures
quoted from the design and research documents checked out** — the same split this project has now
seen on every milestone of both step-4 plans.

**Mi12: the module doc defined a stratum as the three-part key including the read group**, where
the type has two fields and the architecture files fits under `(ReadGroupId, Stratum)`.
**Category:** naming.

### Nits

Missing `#[must_use]` on the free functions and accessors, where the sibling path applies it
consistently; `OFFSET_BUCKETS` named as a plural for a count; "support" as untranslated statistics
jargon; the three scaffold files carrying no `#[cfg(test)]` block where A1's wording asks for one
(deliberate, and now recorded); `bucket_of` reading as a noun; the `ssr`/STR spelling split having
no glossary line.

## 7. Out of scope observations

- **The citation defect behind Open question 1.** `research/…_2026-08-06.md` §6.8 is "What is not
  measured" and there is no §6.4.1; five locations across the spec and architecture cite them for
  figures they do not contain. Pre-existing; the new code comments inherited it verbatim and now
  cite the spec sections that actually hold the numbers instead.
- **`arch/parameter_prepass_ssr.md` cites §4.2 for the monotonicity walk**, which is spec §4.3;
  §4.2 is "How the four numbers are fitted". The shipped code comment cites §4.3 and is right.
- **`src/ng/repeat_catalog/strata.rs` already owns "stratum"** as a raw `(u8, u64)` pair, and its
  module doc says it exists to serve this fit. The two spellings are now permanent — the catalog is
  a peer module and cannot import from a pipeline stage — so every hand-off converts. Recorded in
  `Stratum`'s doc comment; the lift-to-`types.rs` trigger is the milestone where a driver first
  calls `sample_loci_per_stratum`.
- **`examples/` clippy and the 13 rustdoc links**, both pre-existing and both verified so.

## 8. Missing tests to add now

All five landed in the fix pass except the last, which pins a relation no code path uses yet:

1. `ssr_period_accepts_exactly_the_str_scope` — proptest over `0..=1000`, in place of three points.
2. `the_allele_support_is_every_offset_from_the_clip_to_the_limit` — the whole vector at 0, 3, 5, 6
   and 20 reference repeats.
3. `every_i8_offset_maps_into_the_bucket_range_without_going_backwards` — all 256 offsets.
4. `an_offset_renders_with_its_sign_in_both_directions`.
5. `a_strata_map_walks_in_the_catalogs_own_stratum_order` — **deferred to Milestone C**, where the
   conversion between `Stratum` and the catalog's tuple first exists.

## 9. What's good

- The saturation fixture reaches `i8::MIN` and `i8::MAX`, which is exactly where a naive
  `(offset + HALF_RANGE) as u8` overflows — the strongest fixture in the diff.
- `strata_sort_by_period_and_then_by_repeat_count` anti-correlates repeat count with period, so
  both wrong orderings differ from the asserted one; it discriminates on both halves of the key.
- `the_support_survives_a_repeat_count_far_larger_than_an_offset` picks `u32::MAX`, the one value
  that separates a checked conversion from an `as i8` cast.
- Every doc comment carrying a design number carries it correctly — 22 of 22 checked.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
./scripts/dev.sh cargo test --lib --bins --tests --all-features
./scripts/dev.sh cargo test --lib parameter_estimation::ssr
```

Audit trail: `tmp/review_2026-08-11_ng-prepass-ssr-a1a2a3/` (eight per-category files, and the
patch each agent reviewed).
