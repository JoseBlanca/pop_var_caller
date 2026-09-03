# Code Review: ng_psp_mode_c1
**Date:** 2026-09-03
**Reviewer:** rust-code-review skill (orchestrator; six sub-agents in isolated worktrees, eight category checklists between them)
**Scope:** step C1's uncommitted diff — the `generate-psps` subcommand, on f00d56e9
**Status:** Request-changes (all applied — see the fix report)

---

### 1. Scope

- Reviewed: the working-tree diff of plan step C1, exported as `tmp/review_2026-09-03_c1_generate_psps/c1.patch`.
- Against: commit `f00d56e9` + the patch, branch `ng-psp-mode`.
- In-scope: [src/pop_var_caller_exp/generate_psps.rs](../../../src/pop_var_caller_exp/generate_psps.rs) and its [tests](../../../src/pop_var_caller_exp/generate_psps/tests.rs) (both new), plus wiring in `mod.rs`, `cli.rs` and [src/main_exp.rs](../../../src/main_exp.rs).
- Out of scope: `run_ground.rs` and the `ng::run` fixture module (committed earlier in this milestone), `gatherer.rs` (reviewed at B1).
- Categories (8): reliability, errors, defaults, naming, idiomatic, smells, module_structure, extras, refactor_safety — every "always" category plus defaults (a public command surface), module_structure (multi-file), extras (diff-matches-intent).

### 2. Verdict

**Request-changes.** The command works and its shape is right, but the review found two live defects and a set of test holes that mutation testing proved were not theoretical: **41 mutations were run across four harnesses and 21 survived.** The two live defects are a path built from unvalidated file metadata, and a re-walk that destroys the psp it is replacing.

### 3. Execution status

- Every agent detached at `f00d56e9`, applied the patch and verified branch-only files first.
- Mutation totals: **41 run, 21 survived, 19 killed, 1 changed-no-behaviour** (reliability 14/9, refactor_safety 8/2, defaults 9/4+1, errors 5/4, idiomatic 2 verified rewrites). The reliability agent wrote seven probe tests *before* mutating and ran them green on the clean tree, so every survivor it reports is proven to change behaviour.
- Orchestrator-side: `cargo fmt --check` clean, `clippy -D warnings` exit 0, `cargo test --lib 'pop_var_caller_exp'` 100 passed. The command was also run for real on a tomato CRAM slice: `SRS3394712.psp`, 914,715 bytes, 3.0 s, exit 0.

### 4. Open questions and assumptions

1. **What should a sample name that is not a file name do?** Refuse, or encode? Taken as *refuse, at the door* — the walk is minutes and the alternative is finding out at the write.
2. **Should the read filters and locus-generator settings become flags?** Both commands hard-code them today. Left as-is: making them typeable is a change to both surfaces and a design question for the owner, not this step's. Affects M7.

### 5. Top 3 priorities

1. **B1/B2** — the two live defects: a psp path built from `@RG SM` without checking it, and a stopped re-walk truncating the good psp it was replacing.
2. **B3/B4** — the two paths whose breakage produces wrong results silently and which no test covers: the catalog-against-reference digest check, and the mid-cohort failure.
3. **M1/M2** — four flags never given a non-default value by any test, and `--min-purity` missing from the cross-command comparison that exists to stop the two modes drifting.

### 6. Findings

#### Blocker

- `generate_psps.rs:364` — **B1: the psp's file name is built from `@RG SM` with nothing checking it**
- **Categories:** smells, reliability, errors (convergent, three agents)
- **Confidence:** High
- **Problem:** `output_dir.join(format!("{sample}.psp"))` takes the sample name verbatim from the alignment header, which `read_groups.rs` passes through `from_utf8_lossy` and its own doc calls "a *label*". `SM:../elsewhere` writes outside `--output-dir`; an absolute `SM` discards `--output-dir` entirely; `SM:lane/1` fails at the write with the sample's whole walk already decoded; an empty `SM` yields a hidden `.psp`.
- **Fix:** refuse a name that is not exactly one normal path component, for every sample, before the reference is read.

- `generate_psps.rs:307` — **B2: a re-walk that stops destroys the psp it was replacing**
- **Categories:** errors
- **Confidence:** High
- **Problem:** `write_psp` is handed the final path and `PspWriter::create` truncates. The command's advertised repair is *re-run the one sample that failed, into the same `--output-dir`* — so a second failure truncates the good psp at the first byte and leaves a stump `PspReader::open` refuses. The writer's own doc hands the choice to the caller: "write to a new path and rename if that matters".
- **Fix:** walk into `<sample>.psp.partial` and rename once whole; remove the stump on failure.

- `generate_psps.rs:274-278` — **B3: the one check that catches a catalog built on another reference has no test**
- **Categories:** reliability (mutation 16)
- **Confidence:** High
- **Problem:** only the FASTA-verified view of the reference carries digests, and only digests catch a catalog built on a *different* reference with the same contig names and lengths. Passing the `.fai`'s digest-free view instead left all 100 tests green. A psp written against the wrong catalog puts every repeat tract at the wrong coordinates genome-wide, in a file that opens cleanly.
- **Fix:** the agent's verified test — a second reference of the same shape with different bases, a catalog built on it, `--catalog` pointed at that.

- `generate_psps.rs:290-314` — **B4: the mid-cohort failure path is not tested at all**
- **Categories:** reliability (mutations 14, 8, 5), errors, extras
- **Confidence:** High
- **Problem:** three breakages survive: swallowing a stopped walk (**exit 0 with a psp missing**), blanking the sample name and path from the message, and walking samples in sorted rather than first-seen order — the last being a documented invariant whose only observable consequence lives on this path.
- **Fix:** one test with a third sample whose file cannot be opened, asserting the error names it and that the earlier samples' psps are finished.

#### Major

- **M1: four flags are never given a non-default value** (`--regions`, `--catalog`, `--min-purity`, `--build-index-if-missing`) — severing any from the run survives the suite (reliability, mutations 1/2/3/9). The fixture hides the catalog case by writing it exactly where the default looks. A dropped `--regions` writes a whole-genome psp that is internally consistent and silently the wrong file.
- **M2: `--min-purity` is absent from the cross-command defaults comparison** (defaults, refactor_safety — convergent, mutation-proven twice). The tuple has four terms; moving that one default in `generate-psps` alone leaves the suite green, and a drift there makes a later calling run refuse every psp.
- **M3: the provenance timestamp silently falls back to the Unix epoch** (errors, defaults, reliability, idiomatic — four agents). The two other provenance sites in the tree refuse with a `TimestampFormat` error. The timestamp is the field §12.1's byte-identity oracle exempts, so a constant one makes two runs identical for the wrong reason.
- **M4: the output directory is judged last, not first** (smells, extras, errors) — after the reference read, the segmentation and every alignment header, contradicting the function's own comment three lines above and the sibling command's explicit ordering.
- **M5: the `.psp` extension is unpinned** (reliability, mutation 4) — every test locates the file by calling `psp_path_for` again, so dropping the extension leaves them all green.
- **M6: no command-level test walks a sample with any reads** (reliability, defaults) — both fixture BAMs hold zero records, so nothing can tell a psp holding a sample's evidence from a well-formed empty one, and two of the defaults agent's mutations were uncatchable for that reason alone.
- **M7: eleven hard-coded behavioural values are invisible at the command surface** (defaults) — the read filters and the five locus-generator knobs. The filters are recoverable from the psp header; of the generator's five only `max_record_span` is.

#### Minor

Five of six error variants never constructed by a test (reliability); `files_of` re-deriving by string scan the grouping the table already holds (idiomatic, smells, module_structure — convergent); the provenance subcommand literal untied to the clap variant, so a rename leaves every psp naming a dead subcommand with the suite green (naming); `Walk` covering both "no file was created" and "a partial file is on disk", two states needing different recovery (naming); the 22-line reference-open block duplicated verbatim with the sibling command, which `run_ground`'s own doc claims to own (smells); the five routing `#[arg]` blocks byte-identical to the sibling's, which `#[command(flatten)]` over the existing `RepeatRouting` would collapse (smells, Medium confidence); the `a_cohort_on_disk` fixture duplicated 49 lines (smells); **`gatherer.rs:233-235` asserting in the present tense that `generate-psps` checks for an existing psp before truncating, which it does not and C3 has not landed** (three agents, cross-category); `--regions`' help pointing at `call-from-psps`, which does not exist; `--max-str-len` promising a report this command does not print; the `# Errors` paragraph naming three pre-open checks where there are four; "psp" never defined on any CLI surface a geneticist meets; `psp_path_for` `pub` with no caller outside its module; a test message claiming the whole reference is analysed while asserting only non-emptiness; the two-psp equality assertion true by construction; an empty output directory left behind by a failed run; `GeneratePspsArgs` lacking the `Clone` its sibling derives.

#### Nits

`GeneratePspsArgs`' `--min-copies` default a third hand-written copy of the six numbers; four long-lived locals named by shape (`info`, `verify`, `analysed`, `with_checksums`), all inherited from the sibling so worth renaming only in both files at once.

### 7. Out of scope observations

- The `.fai`-vs-verified reference distinction is subtle enough that the sibling command has the same shape and the same absent test; worth a look when `call-from-psps` lands.
- `ground_request` is duplicated verbatim between the two commands — it would disappear with the `#[command(flatten)]` change.

### 8. Missing tests to add now

The reliability agent supplied seven complete bodies, all run green against the unmutated tree (`107 passed`): the catalog-on-another-reference refusal, the mid-cohort stopped walk, `--regions` narrowing the ground, `--catalog` read rather than the sibling file, `--min-purity` reaching the criteria, `--build-index-if-missing` reaching the open, and the psp names read from the directory. Plus: a sample-name refusal, a stopped-re-walk-preserves-the-old-psp test, an unreadable alignment naming itself, a subcommand-name tie, and one that asserts a sample's reads reach its psp.

### 9. What's good

- `files_of`'s deduplication and ordering are genuinely pinned — three separate breakages each killed by a named test, and the dedup test carries a guard assertion against fixture rot.
- Every struct literal in the diff is exhaustive (`SampleWalkInputs` 5/5, `WriterProvenance` 8/8, `GroundRequest` 4/4, `RepeatRouting` 5/5, `GeneratePspsArgs` 11/11 in both test literals): a field added anywhere compile-breaks the command.
- The diff is C1's scope exactly — C2's report and C3's overwrite guard are both absent, and the extra flags beyond C1's list are what `segments_over` requires rather than additions.
- The error enum follows the sibling's shape variant for variant, and the shared ground refusals are carried transparently so one mistake reads identically from either command.

### 10. Commands to re-verify

`scripts/dev.sh cargo fmt --check`; `… clippy --all-targets --all-features -- -D warnings`; `… cargo test --lib 'pop_var_caller_exp'`; and the real run: `generate-psps --reference <ref> --catalog <cat> --alignment <cram> --regions <bed> --output-dir <dir>`.

### Author response convention
Address findings by identifier (B1–B4, M1–M7) in the fix-application report.
