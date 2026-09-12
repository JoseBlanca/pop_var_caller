# Promoting ng to the production caller

*Status: draft, 2026-09-12. Branch `promote-ng`.*

This plan turns a settled decision into build order. It is **not** a place for new design: no
module is redesigned, no behaviour changes, and every step is proven by the caller writing the
same file it wrote before the step.

**The decision (owner, 2026-09-12).** ng has shown it works, so the production caller goes and
ng becomes the one caller in the crate. The owner's three-step proposal — delete production,
check ng, move `src/ng/` to `src/` — is reordered here for one reason: **ng's shipped code
imports production in 77 places**, so deleting production first leaves nothing that builds.
The order that works is *sever, freeze the oracles, delete, move*; §3 says why each stands
alone. The owner accepted that order and the oracle-freezing step in the same conversation.

Design sources this plan cites: the freeze paragraph at the head of
[`src/ng/mod.rs`](../../../src/ng/mod.rs) (what ng may and may not reuse from production),
[`doc/devel/ng/arch/module_layout.md`](../ng/arch/module_layout.md) §"Where ng lives" (ng is a
module of this crate so that the promotion is an in-crate move), and
[`run_streaming.md`](../ng/spec/run_streaming.md) §12.3 (the byte-for-byte VCF comparison
this plan uses as its oracle at every checkpoint).

## 1. Vocabulary

Three words carry the whole plan, so they are fixed here.

- **Production** — the two-phase caller shipped as the `pop_var_caller` binary: `src/main.rs`,
  `src/pop_var_caller/`, and the modules only it runs — `var_calling`, `vcf`, `pileup`, `psp`,
  `paralog`, `sample_summary`, `baq`, `ssr`, `norm_seqs`, `genetics.rs`, `pileup_record.rs`.
  102,218 lines of `src/` outside ng, about 1,520 tests.
- **ng** — `src/ng/` (the caller) and `src/pop_var_caller_exp/` + `src/main_exp.rs` (its
  command surface, the `pop_var_caller_exp` binary). 321,589 lines, 4,332 + 209 tests.
- **Shared infrastructure** — modules that are neither caller: `bam` (alignment-file reading),
  `fasta` (reference reading), `regions`, `error_render`, `iter_ext`. ng is already their
  heaviest consumer (`fasta`: 18 ng files against 12 in `pileup`). **They stay where they are
  and are not touched by this plan.**
- **An oracle test** — a test in `src/ng/` that runs production's code and asserts ng's copy
  agrees with it. Thirteen files, 103 tests, about 11,900 lines. Two of them (`copy_fidelity.rs`
  in `paralog/` and in `locus_generation/pileup/`) compare production's *source text* through
  `include_str!`; the other eleven compare *computed results*.

## 2. Scope

**In**

- ng stops importing anything from production (77 import sites, 8 modules — table in §5.B).
- Every oracle test either becomes a fixture test — production's answers written to a file once,
  the test comparing against the file — or is deleted with a stated reason (§5.C).
- Production is deleted: the eleven modules above, `src/main.rs`, its 7 integration tests, its
  examples (about 30 of 86), its 6 criterion benches, the scripts that drive its binary.
- ng moves to `src/`; `pop_var_caller_exp` becomes `pop_var_caller`; subcommand names are unchanged.
- The two clippy `allow`s `Cargo.toml` carries only for production
  (`chunks_exact_to_as_chunks`, `manual_slice_fill`) go, and ng is made clean under them.

**Out** (each handed to §8)

- Pruning `bam`, `fasta`, `regions` of items only production used.
- Moving `doc/devel/ng/` up, and repairing the several hundred `src/ng/...` paths in Markdown.
- Renaming the `ng_` prefix off examples, tests, benches and scripts.
- Rewriting PROJECT_STATUS's protected *About this project* paragraph, which describes
  production's pipeline (owner's edit).

## 3. Principles that fixed the order

1. **Sever while the oracles are still alive.** Rehoming `PreparedRead`, the CIGAR vocabulary
   and the left-aligner into ng is the one part of this plan that can go quietly wrong — a
   mate role or an alignment end off by one is a wrong genotype, not a panic. Doing it *before*
   deletion means `left_align_parity` and `locus_generation/pileup/parity` still run production
   beside the rehomed code and catch a slip. Deleting first would throw the guards away in the
   same commit that needs them.
2. **Freeze before you delete.** An oracle test's value is the comparison, not the code it
   compares against. Writing production's answers to a fixture keeps the comparison at the
   cost of one pass over thirteen files. Text guards (`include_str!`) are the exception: they
   assert "identical to the original", which has no meaning once the original is gone.
3. **A deletion commit deletes.** Every repoint, rename and doc-link repair lands before the
   deletion, so the deletion diff is `-` lines only and a reviewer can read it as a list.
4. **A rename commit renames.** Moving 291,000 lines from `src/ng/` to `src/` touches every
   file; if that commit also changes a line of behaviour nobody can find it. The move is one
   commit with `git mv`, so history follows the files.
5. **One oracle at every checkpoint: the VCF, byte for byte.** Nothing in this plan is meant to
   change a call. So the same cohort called by the pre-plan binary and by the binary after each
   milestone must write identical files — the check `run_streaming.md` §12.3 already defines and
   `scripts/ng_mode_equivalence_oracle.sh` already runs, with its one exemption: the
   `##commandline` header line, which names the binary and so changes at §5.E by construction.
6. **Container builds**, via `./scripts/dev.sh cargo …`, as everywhere in this repo.

## 4. Preconditions (checked before step A1)

- The `promote-ng` worktree exists at `../pop_var_caller-promote-ng`, branched from `main` at
  `fa44d682` or later. ☑ (created 2026-09-12)
- `./scripts/dev.sh cargo test` is green on that commit, and the suite count is recorded in
  `tmp/promote_ng/baseline/test_count.txt`.
- A release build of `pop_var_caller_exp` from that commit exists, and the inputs for the
  oracle are on disk: the tomato cohort's CRAMs under `benchmarks/tomato2/bams/` (33–98 MB
  each), the reference under `$HOME/genomes`, a repeat catalog and a regions file that
  `ng_mode_equivalence_oracle.sh` accepts. Four CRAMs over one chromosome are enough — the
  oracle is identity, not coverage.

## 5. The steps

Every step compiles and passes `cargo test` before the next starts. Steps marked **own commit**
must not be bundled with a neighbour, so a `git bisect` on the oracle can name them.

### Milestone A — the baseline the rest is measured against

- ✅ **A1. Write the baseline.** Build release; run `call-from-alignments` and
  `generate-psps` + `call-from-psps` over the four-CRAM cohort into
  `tmp/promote_ng/baseline/`; record the md5 of each VCF with `##commandline` stripped, of the
  psps' bodies (header timestamp stripped, as spec §12.1 allows), and of the parameters file
  `estimate-parameters` writes from those psps. Record `cargo test` totals for the library, the
  integration tests and the examples separately.
  *Depends:* preconditions. *Source:* `run_streaming.md` §12.1, §12.3.
- ✅ **A2. Freeze the oracle inputs in a script**, `scripts/promote_ng_oracle.sh`, that takes a
  binary path and writes the same md5s to a directory, so every checkpoint below is one command
  and one `diff`. It wraps `ng_mode_equivalence_oracle.sh` rather than repeating it.
  *Depends:* A1. *Source:* `scripts/ng_mode_equivalence_oracle.sh`.

> **Checkpoint A: the baseline exists and reproduces itself.** Running A2 twice against the same
> binary gives identical md5s (this catches a timestamp or a path leaking into a file). Pause for
> review.

### Milestone B — ng imports nothing from production

The 77 import sites, grouped by what they import. Each row is one step; each step copies the
item into ng, repoints every ng importer, and leaves production's copy untouched — the freeze
rule in `src/ng/mod.rs` still holds until Milestone D. Where ng already has its own version of a
type and converts to production's (`PreparedRead`, `MateRole`), the step deletes the conversion
and makes ng's version the only one.

| step | what ng imports today | sites | where it goes in ng |
|---|---|---:|---|
| B1 | `pop_var_caller::common::format_md5_hex` (a 10-line hex formatter) | 6 | `ng/reference_info.rs`, beside the md5 it formats |
| B2 | `psp::varint`, `psp::errors::VarintError`, `psp::index::checksum_index` | 10 | `ng/psp/varint.rs` (copied whole, 314 lines with its tests) and `ng/psp/index.rs` |
| B3 | `pileup_record::ChainId` | 10 | `ng/types.rs` — it is the id vocabulary the psp, the walker and the merge share. `benches/ng_psp_perf.rs` and `examples/dhat_ng_psp.rs` import it too and are repointed here |
| B4 | `genetics::{lgamma, PROBABILITY_FLOOR, MIN_ALT_CONCENTRATION}` (and `ALPHA_REF`, `alpha_from_diversity` in tests) | 9 | `ng/genetics.rs` — the five items and their tests, nothing else of the 595 lines |
| B5 | `pileup::walker::CigarOp`, `WalkerConfig`, the four `DEFAULT_*` constants | 20 | `WalkerConfig` and the constants to `ng/locus_generation/pileup/mod.rs`, the only consumer. **`CigarOp` went to `src/bam/alignment_input.rs`, not into ng** — see the note below |
| B6 | `pileup::walker::indel_norm::left_align_indels` and the `norm_seqs::normalize_alleles` it calls (464 + 214 lines) | 1 | `ng/alignment/left_align_structured.rs`, which is documented as "1a *is* production's `left_align_indels`" |
| B7 | `pileup::walker::{PreparedRead, MateRole}` and `baq_engine::prepare_passthrough` with its four helpers | 5 | `ng/read/prepared_read.rs` already defines ng's own `PreparedRead` and `MateRole`; `prepare_passthrough` is rewritten to build ng's directly, and `from_production` / `into_production` are deleted. **Own commit** — the oracle is `left_align_parity` (still live) plus Checkpoint B's VCF |
| B8 | the 39 doc-comment links into production (`[`…`](crate::pileup::…)`) | 39 | repointed to ng's copy, or turned into plain text where the link was historical. `cargo doc` with `-D warnings` is the check |

- ✅ **B1** · *Depends:* — · *Source:* `src/ng/mod.rs` freeze paragraph ("copies the code").
- ✅ **B2** · *Depends:* — · *Source:* `psp_file_format.md` (the LEB128 varint the record and index encodings use).
- ✅ **B3** · *Depends:* — · *Source:* `psp_chain_id_encoding.md`.
- ✅ **B4** · *Depends:* — · *Source:* `calling_priors.md` (the constants' meaning).
- ✅ **B5** · *Depends:* — · *Source:* `module_layout.md` (`read/` owns the decoded read).

  **Deviation, recorded 2026-09-12: `CigarOp` is declared in `src/bam/alignment_input.rs`,
  not in ng.** The plan put it in `ng/read/aligned_read.rs`, and building it showed why that
  is the wrong home: the input stage *both* produces a CIGAR run (`cigar_to_ops`) and reads
  it back (`cigar_ref_span`, `cigar_is_bad`), and `MappedRead::cigar` is a field of that
  type. `bam` is shared infrastructure this plan keeps, so declaring the type in ng would
  have left `bam` importing a second one from the tree being deleted, and forced a conversion
  at every crossing — 18 lines of mapping existing only to undo a misplacement. Declared in
  `bam` there is one type, no conversion, and the parity oracles keep typechecking.
  `pileup/walker/mod.rs` re-exports it in one line so production compiles unchanged; that
  changes nothing production computes, and the file goes at D3 regardless.
- ✅ **B6** · *Depends:* B5 (`CigarOp`) · *Source:* `alignment.md` (algorithm 1a, the structured left-aligner).
- ✅ **B7** · **own commit, do not bundle** · *Depends:* B5 · *Source:* `read_preparation.md` §1.
- ☐ **B8** · *Depends:* B1–B7 · *Source:* `Cargo.toml` `[lints.rustdoc]` (`broken_intra_doc_links = "deny"`).
- ☐ **B9. The forward guard.** Add to `scripts/precommit-check.sh`, beside the existing
  "production must not depend on ng" step, its mirror: outside `#[cfg(test)]` modules, `src/ng`
  must not name `crate::{pileup, psp, pileup_record, genetics, pop_var_caller, var_calling, vcf,
  ssr, paralog, sample_summary, baq, norm_seqs}`. The compiler is the real check — production is
  deleted in D — but the guard makes B's end state visible before D starts.
  *Depends:* B8. *Source:* `scripts/precommit-check.sh` step 1.

> **Checkpoint B: ng is self-contained and calls the same.** A2's md5s equal A1's; `cargo test`
> totals equal A1's (nothing added, nothing lost — the copies bring their tests, the conversions
> take theirs). Pause for review.

### Milestone C — the oracle tests keep their comparisons and lose their dependency

One step per file. The rule applied to each: **if the test computes production's answer at test
time from inputs the test builds, freeze that answer into a fixture and compare against the
fixture; if it compares source text, delete it; if it uses production only as a helper, copy the
helper.** Fixtures live beside the test in a `testdata/` directory as the module's existing
fixtures do (`calling/allele_candidates/testdata/`, `calling/parameters_file/testdata/`,
`read/input/testdata/`), are under 1 MB each
— a randomised differential is frozen at a fixed seed and a case count that fits — and each
carries a header line naming the production commit that wrote it.

| file | tests | what it compares | action |
|---|---:|---|---|
| `locus_generation/pileup/parity.rs` | 14 | ng's walker records against production's, same read stream | **freeze** — but first confirm what remains: PROJECT_STATUS (2026-09-11) records that the two *whole-output* differentials were retired because ng's records are deliberately no longer production's. Anything left asserting equality on a record ng has since changed is deleted with that citation |
| `locus_generation/pileup/copy_fidelity.rs` | 1 | source text of 4 walker files | **delete** |
| `paralog/copy_fidelity.rs` | 10 | source text of 4 filter files | **delete** |
| `paralog/production_parity.rs` | 14 | scorer ratios, bit for bit, on randomised loci | **freeze** |
| `calling/genotype_table_parity.rs` | 4 | genotype enumeration order, coefficients, homozygous rows | **freeze** — small, exact, and three of its four failures are silent |
| `calling/loop_parity.rs` | 9 | frequency loop, genotype for genotype, on a shared table | **freeze** |
| `calling/quality_parity.rs` | 16 | site quality at the same prior | **freeze** |
| `calling/allele_candidates/ssr_production_differential.rs` | 5 | production's candidate set on 269 tomato tracts | already fixture-driven (`testdata/tomato_tract_*.csv`); **copy** the three production rules it re-implements in-file, drop the imports |
| `read/left_align_parity.rs` | 5 | `PreparedRead` fields against production's `process_read` | **freeze** |
| `alignment/delimit_parity.rs` | 3 | tract offsets against production's `delimit_read` | **freeze** |
| `alignment/leftmost_property.rs` | 20 | a *definition* (leftmost-ness), not production | **copy helpers** — it only uses `CigarOp` (now ng's after B5) and mentions `normalize_alleles` in a doc link |
| `scanner_parity.rs` | 1 | ng's scanner against a golden catalog already committed under `tests/data/tandem_repeat/` | **copy** production's post-filter `build_loci` and `TrfRecord::for_test` into the test; the golden file is already the fixture |
| `window_coverage/production_parity.rs` | 1 | window means against `SlidingWindowCoverageAccumulator` on generated streams | **freeze** |

- ☐ **C1–C13**, one checkbox per row, in table order. Each: *Depends:* B9 · *Source:* the
  file's own header, which names the spec section it guards.
- ☐ **C14. The four parity examples** — `ng_psp_against_production.rs`, `ng_psp_parity.rs`,
  `ng_psp_head_encoding.rs` (`test = true` in `Cargo.toml`), `paralog_score_parity.rs` — are
  deleted with a line each in the report saying which document already holds their result.
  *Depends:* — · *Source:* PROJECT_STATUS entries for each.

> **Checkpoint C: production is now imported only by production.** `grep -rn 'crate::\(pileup\|psp\|…\)' src/ng` is empty; `cargo test` runs the frozen tests green; A2's md5s equal A1's.
> Pause for review — this is the last checkpoint where the owner can ask for another oracle to be
> frozen, because after D there is nothing to freeze it from.

### Milestone D — delete production

Leaves first, so each commit compiles. `lib.rs` loses its `pub mod` line with each tree.

- ☐ **D1. The command surface and its tests.** `src/main.rs`, `src/pop_var_caller/`, the
  `[[bin]] pop_var_caller` entry; `tests/{cohort_cli,pileup_cli,psp_to_pileup,thread_budget,
  contamination_estimation}_integration.rs`.
  *Depends:* C. *Source:* §2 In.
- ☐ **D2. The calling engine.** `src/var_calling/`, `src/vcf/`;
  `tests/{cohort_vcf_writer,posterior_engine}_integration.rs`; benches
  `cohort_var_calling_perf`, `paralog_scoring_perf`.
  *Depends:* D1.
- ☐ **D3. The per-sample stage.** `src/pileup/`, `src/psp/`, `src/paralog/`,
  `src/sample_summary/`, `src/baq/`; benches `pileup_walker_scaling`, `baq_perf`,
  `psp_reader_perf`, `psp_writer_perf`.
  *Depends:* D2.
- ☐ **D4. The STR caller and the leaves.** `src/ssr/`, `src/norm_seqs.rs`, `src/genetics.rs`,
  `src/pileup_record.rs`.
  *Depends:* D3.
- ☐ **D5. Examples.** The 22 examples that import only production (`dhat_var_calling`,
  `het_baseline`, `psp_rechunk`, `ssr_psp_dump`, …, listed by `grep -l 'pop_var_caller::\(pileup\|var_calling\|…\)' examples/`), and the 8 ng probes that borrowed one production item
  (`ng_depth_term_family` → `lgamma`, `ng_normalizer_screen` → `CigarOp`, …) repointed to ng's
  copies.
  *Depends:* D4.
- ☐ **D6. Scripts and benchmark drivers.** `scripts/{cohort_memory_vs_samples,psp_block_window_sweep}.sh`
  drive the deleted binary — deleted; `scripts/attribute_peak.py` lists production's source
  directories in its heap-attribution table — trimmed. Under `benchmarks/`, the 15 `run_ours_*` /
  `perf_ours_*` / `build_psp` drivers that call `pop_var_caller var-calling` — deleted; their
  *results* and reports stay, since they are the record ng was measured against.
  `benchmarks/lib/common.sh`'s binary discovery moves to the one binary.
  *Depends:* D1.
- ☐ **D7. `Cargo.toml` and the gates.** Remove the two clippy `allow`s and fix what fires in ng
  (`as_chunks::<N>()` for `chunks_exact`, per the comment that asked for it); drop the
  `[[bench]]` entries from D2–D3; `precommit-check.sh` step 1 (production must not import ng)
  and B9's mirror both become vacuous — replace them with nothing.
  *Depends:* D4. *Source:* `Cargo.toml` `[lints.clippy]` comment ("Drop both `allow`s when
  production is retired").
- ☐ **D8. `src/ng/mod.rs`'s header** — the freeze paragraph and the oracle inventory describe a
  world that no longer exists; rewrite to say what ng *is*, with a dated line saying production
  was deleted here. `src/lib.rs`'s crate doc likewise.
  *Depends:* D4.

> **Checkpoint D: one caller, same calls.** `cargo build`, `clippy --all-targets -D warnings`,
> `test`, `doc -D warnings`, `bench --no-run` all green; `cargo test` total equals A1's minus
> production's tests and the deleted oracle tests, and the report says the two numbers;
> A2's md5s equal A1's. Pause for review.

### Milestone E — ng moves to `src/`, and its binary takes the name

- ☐ **E1. The move.** `git mv src/ng/<module> src/<module>` for each of ng's 19 top-level
  entries; `src/ng/mod.rs`'s declarations and re-exports fold into `src/lib.rs`; every
  `crate::ng::` becomes `crate::` and every `pop_var_caller::ng::` in `tests/`, `examples/`,
  `benches/` becomes `pop_var_caller::`. No module name collides — production's `psp`, `vcf`,
  `paralog` are gone at D3–D4, which is why E follows D. **One commit, rename-only.**
  *Depends:* D8. *Source:* `module_layout.md` §"Where ng lives".
- ☐ **E2. The binary.** `src/main_exp.rs` → `src/main.rs`; `src/pop_var_caller_exp/` →
  `src/cli/`; clap `name = "pop_var_caller"`; the `[[bin]]` entry; `##commandline` now says
  `pop_var_caller`. The 4 `scripts/ng_*.sh` and 6 `benchmarks/**/run_ng*` drivers that look
  for `pop_var_caller_exp` in `target*/release/` are repointed.
  *Depends:* E1. *Source:* §2 In ("subcommand names are unchanged").
- ☐ **E3. CI.** `.github/workflows/ci.yml`'s release-test step is scoped to `ng::calling`;
  rescope to `calling`. Nothing else in CI names a path.
  *Depends:* E1.

> **Checkpoint E: the same file under a new name.** A2's md5s equal A1's — `##commandline` is the
> exempted line and the only one that may differ. Pause for review.

### Milestone F — close

- ☐ **F1. PROJECT_STATUS** — the *Last completed task* entry with the numbers from every
  checkpoint; the *About* paragraph flagged for the owner (§8).
- ☐ **F2. CLAUDE.md** — the "Source tree" layout pointer still names `src/main.rs` as
  production's CLI; and the assistant's own memory notes about the production freeze are
  retired.
- ☐ **F3. Implementation report** under `doc/devel/reports/implementations/`, with the
  before/after line and test counts and the md5 table.

## 6. Verification summary

| milestone | proven by |
|---|---|
| A — baseline | A2 run twice against one binary gives identical md5s |
| B — sever | md5s = A1; test totals = A1; `left_align_parity` and `pileup/parity` green with production still present |
| C — freeze | md5s = A1; frozen tests green; `grep` for production paths in `src/ng` outside `#[cfg(test)]` is empty |
| D — delete | md5s = A1; all five cargo gates green; test total = A1 − production − deleted oracles, both numbers stated |
| E — move | md5s = A1 except `##commandline`; `git log --follow` on a moved file reaches its ng history |

## 7. Sizes, so the report can be checked

| | before | expected after |
|---|---:|---:|
| lines under `src/` | 423,807 | about 325,000 |
| library tests | 4,332 (ng) + 209 (exp) + ≈1,520 (production) | ≈4,450 — ng's, minus the 11 text-guard tests, plus fixtures' |
| examples | 86 | about 60 |
| criterion benches | 12 | 6 |
| binaries | 2 | 1 |

## 8. Out of scope — handed on

| item | why not here | where it goes |
|---|---|---|
| Pruning `bam`, `fasta`, `regions` of items only production reached (e.g. the BAQ-related read filters in `alignment_input`) | a `pub` item in a library never fires `dead_code`, so the pruning needs its own inventory; and it changes shared code this plan promised not to touch | a follow-on plan, `shared_infrastructure_prune.md` |
| Moving `doc/devel/ng/{spec,arch,impl_plan,reports}` to `doc/devel/` and repairing `src/ng/...` paths in about 200 Markdown files | doc churn that would bury the code diff; nothing in it is checked by a build | the same follow-on |
| Renaming the `ng_` prefix off 40 examples, 5 integration tests, 5 benches, 8 scripts | cosmetic; better done once the doc paths move with it | the same follow-on |
| PROJECT_STATUS's *About this project* paragraph, which still describes `pileup → .psp → DUST → merger → posterior engine` | protected: "do not edit this paragraph" | the owner |
| Whether a fitted inbreeding coefficient improves calls on GIAB | unrelated; open since 2026-09-11 | PROJECT_STATUS, current focus |
