# Fixes applied: ng_paralog_filter_a2

**Date:** 2026-09-06
**Review:** [ng_paralog_filter_a2_2026-09-06.md](ng_paralog_filter_a2_2026-09-06.md)
**Branch:** `ng-paralog-filter` · reviewed at `a8af3866`, fixes folded into that commit
**Outcome:** 6 Major and 12 Minor; **17 applied, 1 deferred**, all nits applied but one.

## The answer

**The four mutations that survived the differential are now killed by it, each verified with
`production_parity` run alone so `copy_fidelity` cannot mask the result.**

| mutation to ng's copy of the scorer | before | after |
|---|---|---|
| `clamp(0.0, cmax)` → `min(cmax)` — the winsor clamp's lower arm | all five differential tests passed | three fail, `a_negative_copy_number_is_winsorised_at_zero_on_both_sides` among them |
| the `precompute.cohort_size != cohort_size` clause deleted | passed | `tables_built_for_another_cohort_are_neutral_on_both_sides` fails |
| the σ₀ length check weakened from `!=` to `<` | passed | `a_mismatched_sigma_slice_is_neutral_on_both_sides` fails |
| the `configs.is_empty()` clause deleted | passed | `an_empty_carrier_set_is_neutral_on_both_sides` fails |

The last two are spec §6 trap 3 and trap 4 respectively: the first turns a precondition failure
into a plausible score over the wrong number of samples, and the second returns a **`NaN` ratio
from a locus that was scored**, which downstream means *unscored*.

The suite goes from 54 tests to **59**, and the scoped suite from 6,329 passing lib tests to
**6,334**.

## Findings table

| ID | Severity | Decision | Status |
|---|---|---|---|
| M1 — four input classes never drawn | Major | Apply | **Applied** |
| M2 — two rustdoc links into production, one unreachable by the mechanism | Major | Apply | **Applied** |
| M3 — the narrowed-generator guard asserted only `count > 0` | Major | Apply | **Applied** |
| M4 — the stream check rebuilt the stream from hand-typed literals | Major | Apply | **Applied** |
| M5 — the guard's structural branch had no test | Major | Apply | **Applied** |
| M6 — a doc comment claimed a spec trap the scorer cannot reach | Major | Apply | **Applied** |
| Mi1 — `GuardedCopy` does not fit A3's copy | Minor | Defer | **Deferred** |
| Mi2 — `GUARDED` restated `guarded_copies()` | Minor | Apply | **Applied** |
| Mi3 — the divergence message names two causes of three | Minor | Apply | **Applied** |
| Mi4 — eighteen literal spaces; no "the guard is off" | Minor | Apply | **Applied** |
| Mi5 — `AtTheFirstItemDoc` can extract nothing and pass | Minor | Apply | **Applied** |
| Mi6 — `ParalogScore` read field by field | Minor | Apply | **Applied** |
| Mi7 — `src/ng/mod.rs`'s register is stale | Minor | Apply | **Applied** |
| Mi8 — "57 times in 100" omits the absent draw | Minor | Apply | **Applied** |
| Mi9 — a cited precedent is xorshift, not splitmix64 | Minor | Apply | **Applied** |
| Mi10 — the `Repoint` header's "one line" and "at run time" | Minor | Apply | **Applied** |
| Mi11 — `4.0` and `200` written out | Minor | Apply | **Applied** |
| Mi12 — the report's mutation values did not reproduce | Minor | Apply | **Applied** |
| Nits (6) | Nit | 5 Apply, 1 partial | **Applied** |

## What was done, finding by finding

### M1 — the four missing input classes · Applied

The generator now draws copy numbers over `[−1, 9)` rather than `[0, 9)`, so the winsor clamp's
lower arm is compared on every locus. The other three classes cannot be drawn — they are
disagreements between the locus and what is handed alongside it — so `score_through_both` takes a
`HandedToBothScorers` whose σ₀ slice length, table cohort size and carrier set are separate
fields, and three tests set one of them apart from the locus.

**Measured, after:** each mutation above fails `production_parity` alone. The negative-copy-number
mutation now fails three tests rather than one, because the widened stream reaches it too.

### M2 — the rustdoc links, and where repoints may reach · Applied

`locus_score.rs`'s two `` [`SingleCopyCoverageModel`] `` link definitions now name
`crate::ng::paralog::`, so ng's own API documentation resolves to ng's copy. One of the two is a
`//!` line **inside production's module header**, which the guard compares as a byte-verbatim
prefix — so substitutions are now applied to the whole of production's file rather than only to
the content past the header. The three declared repoints are the only lines on which ng's copy and
production's differ: `diff` reports **6 lines, 3 pairs**, all three declared.

**And a repoint may now only turn a path into a path** — asserted, `production_line` containing
`crate::paralog` and `ng_line` containing `crate::ng::paralog` — so the mechanism cannot be used to
sanction a changed constant. That is one of the review's nits, closed for the cases the assertion
covers; the general form stays open for A3.

### M3, M4 — the stream, pinned and derived from one place · Applied

The shape counts are pinned exactly: **`[1598, 2210, 1256, 1841, 6259, 1277]`** of 12,600 drawn
observations — absent, zero-read, alt-exceeds-total, degenerate σ₀, past the winsor cap, below
zero — for the same stated reason the per-size locus counts already were. `count > 0` passes
through a 99-in-100 narrowing of any of these rates; an exact count does not.

Both counts and the sweep now come from one `stream_for(cohort_size)` and from `COHORT_SIZES`,
rather than from a second copy of the seed base and the size.

**The last two pinned values were measured, not predicted.** The first draft of the constant
guessed 6,300 and 1,247 from the draw rates; the test reported 6,259 and 1,277, and those are what
is pinned.

### M5, Mi5 — the guard's own rejections · Applied

`a_repoint_is_applied_where_declared_and_nowhere_else` drives every rejection of the `Repoint`
path on synthetic strings: a declared substitution the copy did not make, a body line edited beside
the sanctioned ones, an **undeclared** path edit, and the structural case — a declared line that
is no longer production's, which now reports `THE GUARD COULD NOT RUN` rather than blaming the
copy. `a_guard_that_can_find_no_copied_content_fails_rather_than_passing` closes Mi5: an
extraction that finds nothing is refused instead of comparing `""` with `""` and printing `ok`.

**Mi4's eighteen spaces were a lost `\` continuation, and finding them took the new test.** The
first attempt to rewrite that message did not apply — `cargo fmt` had already collapsed the
string onto one physical line, so the anchor did not match — and the failure surfaced only when
`a_repoint_is_applied_where_declared_and_nowhere_else` asserted on the marker and got the old text
back.

### M6, Mi8, Mi9, Mi10, Mi11 — the prose · Applied

The all-absent test no longer claims to pin spec §6 trap 4: the scorer's refusal is `0.0`, not
`NaN`, and trap 4's sentinel belongs to the run wiring that arrives at step B1. The empty-carrier
test does now touch trap 4's substance from the other side — without the refusal, the marginal is
a log-sum-exp over an empty set and the ratio comes back `NaN` from a *scored* locus — and says so.

"scores nothing 57 times in 100" becomes `1 − (7/8)(3/7) = 5 times in 8`, the σ₀ rate having been
only half of it. `to_toml.rs`'s generator is named as the xorshift it is.
`DEFAULT_MAX_RELATIVE_COPY_NUMBER` and `LOCI_PER_COHORT_SIZE` replace the written-out `4.0` and
`200`. The guard's header says three lines rather than one, and says a rustdoc link reaches
production at *read* time, not run time.

### Mi6 — the differential destructures both scores · Applied

`assert_bit_identical` destructures `production::ParalogScore` and `ng::ParalogScore` rather than
reading fields, so a sixth value added to either is a compile error here instead of a value the
differential quietly stops comparing. This plan's own step B1 requires the same of the spill
encoder.

### Mi2, Mi7 — the two restated lists · Applied

`every_file_in_the_module_is_guarded_ng_s_own_or_released` reads the guarded names off
`guarded_copies()` instead of restating them in a constant whose only cheap agreement was a length
check. `src/ng/mod.rs` names the score alongside the coverage model, and replaces "a test may read
production as an oracle, and two do" — a count that was wrong by six before this commit — with the
rule and the grep that finds the instances, on the same argument the guard's own release table
made before the directory listing replaced it.

### Mi12 — the report's stale numbers · Applied

The implementation report quoted `26.61449637978953` against `26.614488944807874` from a run made
**before** the σ₀ draw rate was fixed, so the stream had moved under it — the failure the
plan-driven skill names as *prose written before its measurement exists*. Re-measured on the tree
as it now stands, the same 1-part-in-10-million perturbation of the winsor clamp fails at
`cohort size 1, locus 0` with `85.12238217765987` against `85.12233325249034`: a difference of
**4.9 × 10⁻⁵ in a value of 85.1, about 1 part in 1.7 million**. Its "7 parts in 10⁹" was wrong on
its own quoted values too, by two orders of magnitude.

## What was deferred, and why

**Mi1 — `GuardedCopy` does not fit step A3's copy.** A3 puts two originals into one ng file:
`src/paralog/prior.rs` and `calibrate.rs:64-118`, a line range in `src/var_calling/`. The type
holds one `production_source` and a `CopyBegins` whose variants both mean "to the end of the
file", `include_str!` is written relative to `src/paralog/`, and three failure messages name that
directory as a literal. The review's proposal — two ng files, one original each, plus an
`ends_before` field and the original's path in the messages — is sound, and building it now would
be building A3's scaffolding during A2 on a shape A3 has not yet met. **It is A3's first task**,
recorded in `PROJECT_STATUS.md` and named in the plan's A3 step when it is reached.

**The nit on `Repoint`'s general form** — that nothing stops a declaration rewriting a constant —
is closed for path substitutions by the assertion above, and stays open in general. The
differential covers `locus_score.rs`; A3 brings a file it does not cover, which is when the
question becomes real.

## Validation

Run in the dev container, on the fixed tree.

    cargo test --all-features --lib "ng::paralog::"
    → test result: ok. 59 passed; 0 failed; 0 ignored; 6290 filtered out

    cargo test --all-features --lib --bins --tests
    → 6334 lib tests passed, 0 failed, 15 ignored
    cargo test --all-features --tests
    → every integration binary green except one test:
      a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

**The lib and integration halves were run separately, and the reason is worth recording.** The
combined run stalled for fifty minutes inside `cohort_cli_integration`'s
`var_calling_regions_bed_restricts_emitted_positions`, which had completed without a slow warning
in all three earlier runs of the same command. Run alone, that binary finishes in **0.07 s** and
all twelve of its tests pass. The cause was machine contention — several container builds from
other sessions at once, each given 16 GB and 8 CPUs — not this branch: every file this step
changes is a `#[cfg(test)]` module under `src/ng/paralog/` or a doc comment, none of which the
production caller's integration tests can see.

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/paralog/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f`, proved so
before A1 began. The gate — no failure `main` does not already have — holds, and the lib count
moves 6,329 → 6,334.

Every mutation above was applied to a file backed up first, and each restoration was confirmed by
`diff` before the next command ran.
