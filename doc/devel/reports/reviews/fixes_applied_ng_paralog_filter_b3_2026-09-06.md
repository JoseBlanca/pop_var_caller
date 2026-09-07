# Fix Application Report: ng_paralog_filter_b3_2026-09-06.md

**Date:** 2026-09-06
**Source review:** [ng_paralog_filter_b3_2026-09-06.md](ng_paralog_filter_b3_2026-09-06.md)
**Source state reviewed against:** `e3e0b6e5` on `ng-paralog-filter`; fixes folded into that commit
**Execution mode:** non-interactive
**Overall status:** Completed

**One thing folded in besides the fixes.** `e3e0b6e5` had a docs commit on top of it (`8b731d41`,
three lines in `PROJECT_STATUS.md` recording that B3's review had been dispatched and not
collected). Collecting the review made that note false, and its content is superseded by this run's
status update, so the two commits are collapsed into one rather than leaving a commit whose message
contradicts the tree. Nothing else was in it.

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 7
- Minors: 15
- Nits: 6 (grouped)

### Outcome totals
- Applied: 19
- Applied with adaptation: 3
- Deferred: 6
- Disputed: 0
- Already fixed: 0
- Failed validation: 0
- Awaiting user answer: 2 (the review's two open questions, raised at Checkpoint B)

### Validation summary

Run in the dev container from this worktree, on the tree being committed.

- `cargo test --all-features --lib "ng::run::paralog_filter::patch"` → `ok. 23 passed; 0 failed`
- `cargo test --all-features --lib "ng::vcf::writer"` → `ok. 22 passed; 0 failed`
- `cargo test --all-features --lib --bins --tests` → **6,439 lib tests passed, 0 failed, 15
  ignored**; one integration test failed
  (`ng_calling_loop_calls_genotypes::a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`),
  which is `main`'s at `a33ada0f`
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` → 3 errors, all
  `needless_lifetimes` in `cohort_merge/build.rs:820`, `:893` and `serial.rs:67` — **the same three
  as `main`, none new**
- `cargo fmt --check` → dirty on **the same 9 unique files as `main`**, none of them this step's
- `cargo doc --no-deps`, `cargo audit`, `--all-targets` → not run, for B1's and B2's reasons
- Performance check → not applicable: no `Apply` touched code reachable from a bench harness
  (`rewrite_filter_and_info` has no caller, and `place_of` is not benched)

**The lib count moves 6,427 → 6,439**: ten tests added to `patch` (13 → 23) and two to
`vcf::writer` (20 → 22).

### Unresolved high-priority findings

None. Every Blocker and Major is Applied or Applied with adaptation.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | the contract is tested only against lines the encoder does not emit | Apply | **Applied** | `patch/tests.rs` | 4 new tests, 12 encoded shapes |
| M1 | Major | the ordering check reads a `place` nothing ties to the line | Apply | **Applied** | `writer.rs`, `paralog_filter/mod.rs`, `writer/tests.rs` | `From<&SpillEntry>`; pinning test |
| M2 | Major | the padding rule has two copies | Apply | **Applied** | `writer.rs`, `encode.rs` | 22 writer tests pass |
| M3 | Major | `RecordPlace.position` is a bare `u64` where the spill entry carries a `Position` | Apply | **Applied** | `writer.rs`, `vcf/mod.rs`, `writer/tests.rs` | 22 writer tests pass |
| M4 | Major | the error's message, doc, `# Errors` and field doc each describe a different condition | Apply | **Applied with adaptation** | `patch.rs` | doc-only; see §4 |
| M5 | Major | two branches have no test; both mutations survive | Apply | **Applied** | `patch/tests.rs` | 2 new tests |
| M6 | Major | nothing says what the writer's state is after it refuses a line | Apply | **Applied** | `writer/tests.rs` | 1 new test |
| M7 | Major | the error names no record | Apply | **Applied with adaptation** | `patch.rs` | doc-only; see §4 |
| Mi1 | Minor | all five rules the module encodes are VCF grammar | Defer | **Deferred** | — | open question 1 |
| Mi2 | Minor | `MISSING` and `PASS` re-spell values the crate owns | Apply | **Applied** | `patch.rs`, `vcf/mod.rs` | 23 patch tests pass |
| Mi3 | Minor | `RecordPlace` not re-exported | Apply | **Applied** | `vcf/mod.rs` | builds |
| Mi4 | Minor | `check_order` reads fields one at a time | Apply | **Applied** | `writer.rs` | 22 writer tests pass |
| Mi5 | Minor | length-check-then-index | Apply | **Applied with adaptation** | `patch.rs` | array form, not slice; see §4 |
| Mi6 | Minor | the owned return makes D2's deferral cost two changes | Defer | **Deferred** | — | see §5 |
| Mi7 | Minor | `write_filter` / `write_info` do not write | Apply | **Applied** | `patch.rs` | renamed to `append_*_column` |
| Mi8 | Minor | `#[non_exhaustive]` missing on the variant | Apply | **Applied** | `patch.rs` | builds |
| Mi9 | Minor | four items `pub` where `pub(crate)` fits | Defer | **Deferred** | — | see §5 |
| Mi10 | Minor | neither refuses a tab or a newline | Defer (behaviour) / Apply (doc) | **Applied with adaptation** | `patch.rs`, `writer.rs` | doc half only; open question 2 |
| Mi11 | Minor | not idempotent, nothing says once-only | Apply | **Applied** | `patch.rs`, `patch/tests.rs` | doc + pinning test |
| Mi12 | Minor | three tests pin inputs the encoder cannot produce | Apply | **Applied** | `patch/tests.rs` | renamed + encoder assertion |
| Mi13 | Minor | `added_bytes` has no test | Apply | **Applied** | `patch.rs`, `patch/tests.rs` | 1 new test, 36 shapes |
| Mi14 | Minor | no property test | Apply | **Applied** | `patch/tests.rs` | 2 proptest properties |
| Mi15 | Minor | two of the five new writer tests cannot fail on the obvious break | Apply | **Applied** | `writer/tests.rs` | whole body compared |
| Nits | Nit | six, grouped | Apply (4) / Defer (2) | **Applied** / **Deferred** | `patch.rs`, `writer.rs` | see §5 |
| §6a | — | five wrong claims in the diff's own prose | Apply | **Applied** | impl report, `patch.rs`, commit message | re-measured |

## 3. Questions asked and answers

None asked during the run. **Two open questions go to the owner at Checkpoint B** rather than
blocking: whether `patch.rs` moves to `src/ng/vcf/` (review open question 1, finding Mi1), and
whether the patch should refuse a tab or a newline rather than documenting that it does not
(open question 2, finding Mi10). Neither blocks C1, and both are recorded in `PROJECT_STATUS.md`.

## 4. Adaptations — where the fix differed from the review's suggestion

- **M4** — the review offered two ways to make the four statements agree: tighten the guard to ten
  pieces, or say nine everywhere. **Taken the second**, doc-only. Requiring a tenth piece would add
  a check against an input the encoder cannot produce (`VcfRecord::new` refuses a record with no
  sample columns), and the guard's job is the pieces the patch must index rather than whether the
  line is a valid cohort record. The variant doc now says so explicitly, so the next reader does
  not re-open it. The uninterpolated `{PIECES}` in the field doc is gone and the bound is stated
  correctly as eight.
- **M7** — the review suggested either adding context fields to the error or making C4 wrap it.
  **Taken the second**, because `rewrite_filter_and_info` genuinely does not know where the line
  came from and inventing a field it cannot fill would be worse than saying whose job it is. The
  variant's doc now names C4 and says what to wrap it with.
- **Mi5** — the review's own measurement showed the two salvaged probes' **slice** pattern buys no
  compiler enforcement (at `PIECES = 10` it compiles and fails 10 of 13 tests at run time) while a
  fixed-size **array** via `try_into` gives `error[E0527]`. **Took the array form.** Verified: the
  code compiles at `PIECES = 9` with all 23 tests passing.
- **Mi10** — split. The doc half is applied on both `rewrite_filter_and_info` and `write_line`,
  naming exactly what is not refused and why. The behaviour half — new error variants for an
  embedded tab or newline — is open question 2, because it puts a scan of every line on a
  per-record path to defend against an input no reachable caller produces.

## 5. Deferred, with reasons

- **Mi1 — moving `patch.rs` to `src/ng/vcf/`.** The review's argument is strong: all five rules the
  module encodes are VCF grammar, and it now imports two of them from `src/ng/vcf`. But deviation 3
  of the implementation report is a recorded decision of a committed step, and reversing one is not
  the implementation loop's call. **Raised at Checkpoint B.** The minimum — importing the
  vocabulary instead of re-spelling it (Mi2) — is applied, so the file no longer holds a second
  copy of what `src/ng/vcf/mod.rs` owns.
- **Mi6 — a `rewrite_filter_and_info_into(&mut Vec<u8>, …)` form.** Real, and the argument that it
  makes D2's deferral cost one change instead of two is correct. Deferred because it adds a second
  public entry point with no caller at all today, to serve a measurement that has not been made.
  Recorded in the implementation report's follow-ups so D2 finds it.
- **Mi9 — narrowing `pub` to `pub(crate)`.** The review's own counterpoint decides it: `spill`,
  `spill_file` and every module on the path are `pub`, so narrowing `patch` alone makes it the odd
  one out. If the surface is worth narrowing it is one decision over the `paralog_filter` subtree
  at the end of the milestone.
- **Two of the six nits** — collapsing `PIECES`'s defining constants (done as a side effect of Mi5,
  which removed both) and `RecordPlace`'s missing `Ord` derives (obtained from `GenomePosition` by
  M3, so nothing to add). Both dissolved rather than being deferred; recorded here so the count
  reconciles.
- **The `VcfWriteError` `#[non_exhaustive]` inconsistency** (a nit, and out of scope) — pre-existing
  on a type this step does not change. Left for whoever touches that error next.

## 6. What the numbers said, re-measured

The review's §6a found five wrong claims about this step's own reach. All five are corrected in the
implementation report, the commit message and the affected doc comments:

| claim | was | is |
|---|---|---|
| how many of the nine failing writer tests pre-date the step | four | **seven** — six calling `write_record`, one calling `write_stream` |
| the patch block's `filtered out` count | 6424 | **6429** (the report's two quoted blocks disagreed by 5) |
| `RecordPlace` versus `SpillEntry`'s head fields | "exactly the three" | two of three matched; **fixed rather than reworded** — `RecordPlace` now carries a `GenomePosition`, so it is exactly the three |
| whether ng's header declares `hiddenParalog` | "ng's header declares" | **it does not**; `FILTER_DECLARATIONS` has five ids and C4 adds the sixth. The doc now says "will declare once step C4 adds" |
| the sample columns as a share of a cohort's line | "all but eight" | **all but nine** — the eight fixed columns plus `FORMAT` |

**And one mechanism, settled before the fan-out and corrected in four places**: ng's encoder cannot
write an empty `INFO` column, because `info_column` pushes `AN=` and `DP=` before any conditional
field. The thinnest `INFO` over every shape ng emits is `AN=0;DP=0` — which
`the_encoder_writes_neither_a_missing_nor_an_empty_filter_or_info` now asserts, so the claim is a
test rather than a sentence.

## 7. Out of scope, recorded not fixed

**`SpilledSample`'s two read-count fields document different rules**
([spill.rs:161-166](../../../../src/ng/run/paralog_filter/spill.rs#L161-L166)). `alt_reads` says
"`0` on every other record"; `ref_reads` says "`AD[0]`" with no qualifier. Spec §3.2 requires
`alt_reads = 0` **and `total_reads = 0`** on every non-biallelic-SNP record, and total is
ref + alt — so a tract spilling its real `AD[0]` would give step C1 a non-zero total and a live
allele term, which is spec §6 trap 1 arriving by a different door. Nothing fills a spill yet (pass
one is C2), so it is latent rather than a bug. It is B1's field, and it bears directly on the
`SpilledSample` shape decision already waiting at Checkpoint B — noted there.

## 8. Commands to re-verify

    ./scripts/dev.sh cargo test --all-features --lib "ng::run::paralog_filter::patch"
    ./scripts/dev.sh cargo test --all-features --lib "ng::vcf::writer"
    ./scripts/dev.sh cargo test --all-features --lib --bins --tests
    ./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
    ./scripts/dev.sh cargo fmt --check
