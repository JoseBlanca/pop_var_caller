# Code Review: ng_paralog_filter_a1
**Date:** 2026-09-06
**Reviewer:** rust-code-review skill (orchestrator) over four category sub-agents, each in its own worktree
**Scope:** commit `1823bba1` on `ng-paralog-filter` — the filter's model constants and per-sample coverage model copied from production into `src/ng/paralog/`, with a textual fidelity guard
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** one commit's diff.
- **Reviewed against:** `1823bba1`, branch `ng-paralog-filter`, parent `a33ada0f`.
- **In-scope files:**
  - [src/ng/paralog/copy_fidelity.rs](../../../../src/ng/paralog/copy_fidelity.rs) — ng's own; the review's centre of gravity
  - [src/ng/paralog/mod.rs](../../../../src/ng/paralog/mod.rs) — ng's own header, declarations and re-exports, plus production's constant block
  - [src/ng/paralog/coverage_model.rs](../../../../src/ng/paralog/coverage_model.rs) — a verbatim copy; its placement, its `use` lines and ng's appended header note only
  - [src/ng/mod.rs](../../../../src/ng/mod.rs) — one module declaration and a header clause
- **Deliberately out of scope:** the *logic* of `coverage_model.rs` and of the copied constant block. They are production's byte for byte and the port's whole claim is that they are unchanged, so a finding asking for an edit there would be a finding against the commit's purpose. Sub-agents were told to route such observations to `Cross-category` as informational. Everything under `src/paralog/`, `src/var_calling/`, `src/sample_summary/` is frozen production.
- **Categories dispatched:** `reliability` (always; and the guard is the only safety net here), `module_structure` (the copy imports a production type — the one structural question the commit raises), `idiomatic` (the guard is the only new Rust), `naming` + `smells` in one agent (the module header and the guard are prose a human reads to decide whether to trust the port).
- **Categories skipped, with reason:** `errors` and `defaults` — every error type and every `DEFAULT_*` in scope is production's, frozen, and un-editable here; findings would be unactionable. `unsafe_concurrency` — no `unsafe`, no threads, no shared state. `tooling` — no `Cargo.toml` change. `extras` — no parser, no untrusted input, no hot path; its "diff matches stated intent" item was run by the orchestrator as step 8a instead. `refactor_safety` — folded into `smells`, which the same agent ran.

### 2. Verdict

**Approve-with-changes.** The commit's central claim — that the copy is production's — held under every check: two independent `diff`s, and five of the eleven mutations run. What the review found is that the *guard* covered less than its prose said, in three ways that each let a real divergence through, and that one forward-looking claim in the header was wrong. All four are fixed; see the fix-application report.

### 3. Execution status

Run by the orchestrator in the dev container, quoted verbatim in each dispatch:

| command | result |
|---|---|
| `cargo test --all-features --lib "ng::paralog::"` | `ok. 28 passed; 0 failed` |
| `cargo test --all-features --lib --bins --tests` | 6,303 lib tests passed; **one integration test failed**, pre-existing |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 3 errors, all `needless_lifetimes` in `cohort_merge/`, pre-existing |
| `cargo fmt --check` | dirty on 9 files, none of them this commit's, pre-existing |

**Commands not run, and why.** `cargo test --all-targets --all-features` does not compile on `main`: `examples/ng_candidate_selection_probe.rs` reads `locus.kind` and hands a `ClosedLocusRanges` to `CohortObservation::over`. It exits 101 having run no test at all, so it measures nothing about any branch. `cargo doc --no-deps` and `cargo audit` were not run; `cargo doc` is already red on `main` with 11 unresolved intra-doc links (`PROJECT_STATUS.md`, *Standing project-wide items*).

**The four pre-existing failures were proved to be `main`'s, not the commit's**, by removing `pub mod paralog;` from `src/ng/mod.rs` — so that none of the commit's files compile — and re-running: 6,275 lib tests pass and the same single integration test fails.

Findings labelled "Needs verification": **0**. Every finding filed below was reproduced by its agent in its own worktree, and the three the orchestrator applied were reproduced again in the working tree before and after the fix.

### 4. Open questions and assumptions

1. **Which fields does ng's `CoverageByGcHistogram` carry?** `window_coverage.md` §7 says the fields the fit does not read are dropped; production's transcribed test fixture builds all eight. Both cannot hold. Affects **M1**. **This is the window-coverage spec's decision, not this plan's** — raised at Checkpoint A, not resolved here.
2. **What makes the temporary import temporary?** Affects **Mi3**. The plan-driven skill forbids this skill from adding a step to the plan, so it is recorded in `PROJECT_STATUS.md`'s open items and raised at Checkpoint A.

### 5. Top 3 priorities

1. **M2** — the guard could go silent: emptying its file list left the whole suite green while a copy had drifted.
2. **M1** — the header tells the next coder the window-coverage swap is one line; it is four, and three of them are inside a file the guard forbids editing.
3. **M3** — 236 copied lines sat in the one file the guard could not reach; deleting both `is_finite()` guards from `GridSpec::new` passed all 6,303 lib tests.

### 6. Findings

#### Major

**M1: [src/ng/paralog/coverage_model.rs:39-45](../../../../src/ng/paralog/coverage_model.rs#L39) — "the swap is one `use` line" is wrong, and the two specs disagree about the histogram's fields**
**Categories:** module_structure (filed), naming (convergent, as a doc-accuracy finding) · **Confidence:** High
The header says the fit reads five of the histogram's fields and ng's copy "drops only fields this file never touches", so swapping the import is the whole change. The fit does read exactly those five. But the file does not end at the fit: its transcribed fixture at [`:714-723`](../../../../src/ng/paralog/coverage_model.rs#L714) constructs all eight, including `n_positions`, `n_skipped_tiles` and `callable_positions`. The reviewer built ng's type as `window_coverage.md` §7 prescribes and got `error[E0560]: struct window_coverage::CoverageByGcHistogram has no field named callable_positions` at `:721`; with the field kept, the swap is one line and the guard fires on exactly it. **Why it matters:** this header is the contract the window-coverage step is planned against, and it promises a one-line commit. Followed as written, that coder hits a compile error inside a file the guard forbids editing.

**M2: [src/ng/paralog/copy_fidelity.rs:52](../../../../src/ng/paralog/copy_fidelity.rs#L52) — nothing checked that the guard's file list still covered every copy**
**Category:** reliability · **Confidence:** High
The list is hand-maintained and the documented release protocol *deletes* from it, so a dropped entry and a deliberate release are indistinguishable. Emptying it while simultaneously drifting a doc line in ng's copy left `cargo test --all-features --lib` at `6303 passed; 0 failed`, with the guard itself reporting `ok` after comparing nothing. The released-files table in the module header was prose nothing read.

**M3: [src/ng/paralog/mod.rs:42](../../../../src/ng/paralog/mod.rs#L42) — 236 copied lines sat in the one file the guard exempted**
**Categories:** reliability (filed), module_structure (convergent) · **Confidence:** High
The header said "This file is ng's own … so it is not guarded", true of the declarations and false of the 236 lines below them, which are production's `src/paralog/mod.rs` byte for byte. The exemption was drawn at file granularity where the divergence is at region granularity. Deleting both `is_finite()` guards from `GridSpec::new` passed the entire `ng::paralog::` suite including `copy_fidelity` — and the four transcribed tests cannot catch it, because they pin values rather than code, which is the exact weakness the guard exists to remove.

**M4: [src/ng/paralog/copy_fidelity.rs:63](../../../../src/ng/paralog/copy_fidelity.rs#L63) — the guard's own comparison had no test**
**Category:** reliability · **Confidence:** High
Its single fixture was the real file pair, currently identical, so every rejecting branch ran only in its accepting direction. Deleting the header-prefix clause and rewriting a production header line in ng's copy together left the whole suite green. A tidy pass that "simplified" the assertion would have converted the project's only cross-tree check into `assert!(true)` with 28 tests still passing.

#### Minor

**Mi1: [src/ng/paralog/mod.rs:8-9](../../../../src/ng/paralog/mod.rs#L8) — the module the header sends the reader to does not exist.** "the run's wiring lives in `crate::ng::run::paralog_filter`", present tense; `grep -rn paralog src/ng/run/` returns nothing. The plan creates it at step B1. **Category:** naming.

**Mi2: [src/ng/paralog/mod.rs:13-15](../../../../src/ng/paralog/mod.rs#L13) — "every file beside this one is production's source byte for byte" is false of one of the two.** `copy_fidelity.rs` is ng's own and says so on its own line 3. This is the paragraph a reader uses to decide which files they may edit. **Category:** naming.

**Mi3: [src/ng/paralog/coverage_model.rs:49](../../../../src/ng/paralog/coverage_model.rs#L49) — the temporary import has nothing making it temporary.** `CoverageByGcHistogram` is a `#[non_exhaustive]` serde type from a TOML document already at version 4 after three breaking bumps, and it appears in the `pub` signature of ng's fit without being re-exported, so a caller must reach into `crate::sample_summary` to name the argument. The guard makes the swap *visible* when it happens; nothing makes it happen. **Category:** module_structure.

**Mi4: [src/ng/paralog/coverage_model.rs:41-42](../../../../src/ng/paralog/coverage_model.rs#L41) — the note cites a spec file that is not in the tree.** `doc/devel/ng/spec/window_coverage.md` exists only on the unmerged `ng-window-coverage` branch. Every other document these headers cite is present, so a reader has no reason to expect this one to be the exception — and this note exists to justify the single most questionable line in the port. **Category:** naming.

**Mi5: [src/ng/paralog/copy_fidelity.rs:52-58](../../../../src/ng/paralog/copy_fidelity.rs#L52) — a 3-tuple whose two same-typed file contents are positionally adjacent.** Transposed, it still passes — line equality is symmetric — and the mistake surfaces only on the day the guard fires, with every message naming the wrong side. **Categories:** smells (filed), naming (convergent: the binding is called `pairs` and holds triples, and that name escapes into both the module header and the failure message as the thing a future author is told to edit).

**Mi6: [src/ng/paralog/copy_fidelity.rs:66,79,88](../../../../src/ng/paralog/copy_fidelity.rs#L66) — two of the three failure messages name no location, and the third misdirects.** The header assert prints no line, so the reader hand-diffs 28 header lines. The length assert runs *before* the per-line loop, so a one-line deletion is reported as a bare count `1157` vs `1158`. The per-line assert fires correctly when *production's* file was the one edited, but tells the reader to release ng's copy, and its `left`/`right` are unlabelled. **Category:** reliability.

**Mi7: [src/ng/paralog/coverage_model.rs:32](../../../../src/ng/paralog/coverage_model.rs#L32) — "byte for byte" promised, line-for-line checked.** `str::lines()` drops a trailing `\r` and cannot see a missing final newline. Stripping ng's final newline (51,397 → 51,396 bytes) and giving one line a CRLF ending both passed the guard. Neither can change what the program computes, so the risk is to the claim — which is what this commit is for. **Categories:** reliability (filed), naming (convergent).

**Mi8: [src/ng/paralog/mod.rs:275](../../../../src/ng/paralog/mod.rs#L275) — the `// non-finite` case does not exercise the finiteness guard.** `GridSpec::new(f64::NAN, 0.6, 40)` is refused by the `lo < hi` clause whether or not `is_finite` is present, because `NAN < 0.6` is `false`. Deleting both `is_finite()` calls left `grid_spec_new_validates_bounds` passing. An *infinite* endpoint is the discriminating input: `GridSpec::new(0.004, f64::INFINITY, 40)` returns `None` on the clean tree and `Some` on the mutant. **Category:** reliability. The test is production's and transcribed, so it may not be edited here.

#### Nits

- The opening paragraph of `mod.rs` stacks four undefined statistical terms in one sentence, in the paragraph that introduces the module to a geneticist, immediately after a sentence that is a model of the house style; and it describes all four as present when only one has landed (`clear-technical-writing` Rule 3).
- `copy_fidelity.rs`'s header asserts the copy is verbatim four times without ever giving its size or its subject, which reads very differently for one file of 1,158 lines than for the sibling guard's eight files of 5,500.
- Five identifier slips in the guard, all the same kind — a name that says the shape or the pronoun rather than the value: `header_lines` returns a count; `their_head`/`our_head` abbreviate "module-header line count"; ng's side is `ours` in one place and `mine` in another; `ours_line`/`their_line` mix plural and singular possessive; `number` is a zero-based index used only as `number + 1`.
- `mod copy_fidelity;` carries no marker saying the file is ng's own, where all three ng-owned test modules in the sibling directory do; and the `pub use` block has no comment saying it reproduces production's surface name for name so that later steps' call sites resolve unchanged.
- The explicit `[(&str, &str, &str); 1]` annotation is unneeded and has to be hand-edited on every add or release — the one edit this file most wants frictionless; and the header-length closure is rebuilt on each loop iteration.
- `grid_specs_have_expected_defaults` overlaps `default_params_match_prototype_constants` on every assertion but one — it reaches `SfsPriorSpec::default()` directly — and would read as duplication to a later tidy.

### 7. Out of scope observations

- **`src/ng/mod.rs:31-36` still says "nothing shipped depends on production."** It was already strained before this commit — shipped ng code imports `crate::genetics` at 7 sites, `crate::pileup::walker::CigarOp`, `crate::pileup_record::ChainId` — and `coverage_model.rs:49` is one more instance. Follow-up: correct that paragraph in whichever branch next touches it.
- **`GridSpec` documents an invariant, enforces it in `new() -> Option<Self>`, and leaves all three fields `pub`**, so a struct literal bypasses it. The doc names this as a deliberate trade. Production's backlog, not this commit's.
- **`window_coverage.md` §7 names `heterozygosity` as a field of `CoverageByGcHistogram` to be dropped.** It belongs to the sibling `HetCounts` (`src/sample_summary/mod.rs:168-180`). That spec's accuracy, on another branch.
- **Running `cargo fmt` on this branch reformats nine files outside this commit.** A contributor who formats before committing picks up an unrelated nine-file diff. Pre-existing on `main`.

### 8. Missing tests to add now

Grouped by what they guard. All were applied; see the fix-application report.

- `the_comparison_accepts_an_appended_note_and_rejects_every_other_edit` — synthetic file pairs driving each rejecting branch of the comparison: a rewritten production header line, a note prepended rather than appended, a dropped header line, a drifted body line, a deleted one, an added one, a stripped final newline, a CRLF ending on a body line and on a header line. Catches M4 and Mi7.
- `every_file_in_the_module_is_guarded_ng_s_own_or_released` — enumerates the module directory and requires every `.rs` file to be classified, and requires the guarded list and the file list to name the same set. Catches M2.
- `the_copies_are_still_productions`, extended to `model_params.rs` — catches M3; verified to fire on the exact mutation that previously passed.
- `grid_spec_new_refuses_an_infinite_endpoint`, as ng's own test beside the transcribed suite — catches Mi8 without editing production's copy.

### 9. What's good

- **The commit does not trust its own tests to prove fidelity, and says why.** `copy_fidelity.rs:9-14` states the reason a green suite cannot establish the property — both trees carry the same transcribed tests — which is the argument that makes the whole guard worth having.
- **The sibling guard's escape hatches were dropped rather than carried.** `is_ng_addition` and `REPOINTS` exist in `locus_generation/pileup/copy_fidelity.rs` for divergences these files do not have; copying them would have been the wrong kind of consistency, and the file says so at the point it omits them.
- **`include_str!` rather than a runtime read**, so a moved or deleted file is a build error instead of a silently skipped case.
- **The `use` line that will change is inside the guarded region**, so the arrangement's own end is forced into a commit that has to explain itself.

### 10. Commands to re-verify

```
scripts/dev.sh cargo test --all-features --lib "ng::paralog::"
scripts/dev.sh cargo test --all-features --lib --bins --tests
scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
scripts/dev.sh cargo fmt --check
```

The last two are expected to fail on `main`'s three lints and nine files respectively; the gate is that neither set grows. `--all-targets` cannot be used until `examples/ng_candidate_selection_probe.rs` compiles again.

Per-category audit trail: `tmp/review_2026-09-06_ng_paralog_a1/{reliability,module_structure,idiomatic,naming_and_smells}.md`.
