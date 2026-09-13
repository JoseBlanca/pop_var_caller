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
  agrees with it. Thirteen files are named for it — 103 tests, 11,927 lines — and two of those
  (`copy_fidelity.rs` in `paralog/` and in `locus_generation/pileup/`) compare production's
  *source text* through `include_str!` while the other eleven compare *computed results*.
  **A sweep before C1 found 46 more such tests in 18 ordinary modules**, plus seven files whose
  dependency on production is not an oracle at all; §5.C's second table lists them and says why
  the count moved.

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
- ✅ **B8** · *Depends:* B1–B7 · *Source:* `Cargo.toml` `[lints.rustdoc]` (`broken_intra_doc_links = "deny"`).
- ✅ **B9. The forward guard.** Add to `scripts/precommit-check.sh`, beside the existing
  "production must not depend on ng" step, its mirror: outside `#[cfg(test)]` modules, `src/ng`
  must not name `crate::{pileup, psp, pileup_record, genetics, pop_var_caller, var_calling, vcf,
  ssr, paralog, sample_summary, baq, norm_seqs}`. The compiler is the real check — production is
  deleted in D — but the guard makes B's end state visible before D starts.
  *Depends:* B8. *Source:* `scripts/precommit-check.sh` step 1.

> **Checkpoint B: ng is self-contained and calls the same.** A2's md5s equal A1's; `cargo test`
> totals equal A1's (nothing added, nothing lost — the copies bring their tests, the conversions
> take theirs). Pause for review.

### Milestone C — the oracle tests keep their comparisons and lose their dependency

**Twenty-six steps, one per file or per group of files that share one production item** — the
first table's thirteen, then the twelve the sweep added. The rule applied to each: **if the test computes production's answer at test
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
| `locus_generation/pileup/copy_fidelity.rs` | 1 | source text of the walker files not yet released — one, `decompose.rs`, by `d9e7b076` | **delete** |
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

- **C1–C13**, one checkbox per row, in table order. Each: *Depends:* B9 · *Source:* the
  file's own header, which names the spec section it guards.
  - ✅ **C1** `locus_generation/pileup/parity.rs` — **14 tests → 8, and one of the 8 frozen.**
    What building found, recorded 2026-09-13: of the 14, three were already `#[ignore]`d under
    the 2026-09-11 retirement (the two whole-output differentials and the real-data one), and
    three more tested only the harness that compared two walkers (the census, the projection,
    the shared reference bytes). Those six were deleted with the citation, as the row allowed,
    and the harness went with them — the file is 1,900 lines, from 4,669. Of the eight that
    stay, **one** still compared ng against production — the malformed-input error stream — and
    it is frozen: production's five outcomes measured at `d9e7b076` and written as literals in
    the test, and a mutation of one literal shown to fail it. Four never touched production and
    three read production's walk only to count what the case generator reaches; those three
    now read ng's walk. Two ng-only checks the retired census ran — at least half as many chain
    ids as reads, and partial-witness runs inside their locus — moved into
    `every_emitted_observation_carries_a_read` rather than going unrun.
  - ✅ **C2** `locus_generation/pileup/copy_fidelity.rs` — deleted. The row said four walker
    files; at `d9e7b076` it guarded one, `decompose.rs`, the others having been released as ng
    changed them. Doc comments in six files and the module's arch doc that described it as a
    running check now say it was retired.
  - ✅ **C3** `paralog/copy_fidelity.rs` — deleted. It held **nine** tests, not the table's ten
    (the tenth `#[test]` the sweep counted is inside a doc comment). Two read production's text —
    one comparing four copies against their originals, one checking the item a span copy stops
    before still exists — and the other seven tested the guard's own machinery (the
    appended-note rule, the repoints, the span markers, the directory bookkeeping). None checked what ng computes: that
    is `production_parity.rs`, frozen at C4. Headers of eight paralog files and `src/ng/mod.rs`
    that described the guard as running are put in the past tense, and `copy_fidelity\.rs`
    leaves `precommit-check.sh`'s exemption list.
  - ✅ **C4** `paralog/production_parity.rs` — **frozen, all 14 tests kept.** Every value
    production's filter returned on the file's generated inputs — the scores of 800 swept
    loci, 200 one-sample loci and 19 hand-built cases; the prior, curve and cut over eight ratio
    streams; the verdict at 320 combinations of stream, target and ratio; the empty and
    unconverged EM; and the fallback at 45 combinations of configuration, histogram and target —
    was recorded at `d9e7b076` by an instrumented copy
    of the test: 6,494 keyed values in `paralog/testdata/production_parity_answers.tsv`
    (274 KB), f64s as bit patterns. The tests still generate the inputs from the same seeds
    and compare ng against the file; each asserts it read every answer in its section. Changing
    one fixture value fails one test; changing one generator draw rate fails three.
  - ✅ **C5** `calling/genotype_table_parity.rs` — **frozen, all 4 tests kept.** Production's
    `shape_for` tables for the 76 shapes the tests compare (27,384 genotype rows) were written at
    `d9e7b076` to `calling/testdata/genotype_tables_production.txt` (733 KB): each genotype as
    its alleles, its log coefficient as a bit pattern, its homozygous allele. The four
    comparisons — count, allele counts in order, coefficient bits, homozygous lookup — read
    from the file. Corrupting one coefficient fails one test.
  - ✅ **C6** `calling/loop_parity.rs` — **frozen, all 9 tests kept.** Production's
    `run_em_columnar` was handed nine likelihood tables across the tests; what it returned on each
    — the genotype called for every sample, its iteration count, whether it converged — is now a
    nine-entry table in the file, keyed by an FNV-1a digest of everything production was handed,
    so a changed input finds no answer and fails. Production's convergence threshold (1e-3), pass
    cap (50) and reference concentration (1.0) are frozen literals ng's constants are checked
    against. Changing one frozen genotype fails one test.
  - ✅ **C7** `calling/quality_parity.rs` — **frozen; 16 tests → 14.** What production's
    `run_em_columnar` returned as site quality on eight distinct likelihood tables, and what
    `vcf::qual_refine::refine_qual` returned on eight sets of reads, are two tables in the file,
    keyed by FNV-1a digests of the inputs; production's shipped pseudocounts `(10, 0.01)` are
    frozen literals. Two tests had no ng half left once production was not run — the
    environment-variable guard on production's strand ramp, and a locus with no alternative
    reads, where ng's correction is never called and the test only read production's answer
    back — and were deleted; a third's production-only assertion now asserts ng's floor
    instead. **Found, not fixed:** `a_homozygous_variant_cohort_is_skipped_by_both` is
    documented as a cohort the correction leaves alone, yet production took its 900 Phred to
    0 and ng agrees within 0.001 — a wrong doc comment or a defect both share; out of this
    plan's scope and raised with the owner.
  - ✅ **C8** `calling/allele_candidates/ssr_production_differential.rs` — **frozen, not copied;
    all 5 tests kept.** *Deviation from the row, recorded 2026-09-13:* the three rules the row
    says to copy were already re-implemented in the file; what still ran production was its own
    selector (`build_rungs` then `assemble_candidates`), used as the answer the re-implementation
    must reproduce. Copying that selector would have meant porting two production modules into a
    test; its answers were recorded instead. What production returned on every tract the tests
    build — 269 tomato tracts with the whole panel, the same with one accession, two hand-built;
    539 distinct, since one tract has a single accession — is in
    `testdata/tomato_tract_production_candidates.csv` (36 KB), keyed by a digest of the whole
    tract. Dropping one candidate from a row fails two tests.
  - ✅ **C9** `read/left_align_parity.rs` — **frozen, all 5 tests kept.** Production's
    `process_read` output (its `--no-baq` arm, `F1` off) for the 8 fixture reads on the uppercase
    and the soft-masked reference — 16 prepared reads, twelve fields each — was recorded at
    `d9e7b076` into `read/testdata/left_align_production_prepared.tsv` and is compared field by
    field. Restoring one row's unshifted CIGAR fails one test. The masked-reference test's check
    that production left every read's CIGAR as the mapper wrote it now checks the recording
    against the fixture's own input, not production's behaviour.
  - ✅ **C10** `alignment/delimit_parity.rs` — **frozen, all 3 tests kept.** Production's
    `delimit_read` answer — a measured byte range, or "ran off an end" — on the 11,986 generated
    cases with a read and the 6 named ones was recorded at `d9e7b076` into
    `alignment/testdata/delimit_parity_production.tsv` (274 KB), keyed by a digest of the case;
    `MAX_SLIP` (10) is a frozen literal. Changing one row's range fails a test. **Lost, and raised
    at Checkpoint C:** the soak that `PVC_PARITY_CASES` drove past 3,000 cases a seed, which is
    what validated `BAND_HEADROOM` (a band one cell too narrow first diverged near case 28,307).
    The docs now say the margin must not be lowered until a production-free reference exists — ng's
    aligner with the band opened to the whole matrix is the candidate.
  - ✅ **C11** `alignment/leftmost_property.rs` — **one doc link repointed; nothing frozen, nothing
    deleted.** As the row said, the file uses no production code: `CigarOp` has lived in
    `bam::alignment_input` since B5, and `normalize_alleles` is named through ng's copy. The one
    production path left was a link to `pileup::walker`'s `indel_norm.rs` for the mismatch-count
    assertion, which ng's own `alignment::indel_norm` carries; the link now points there.
  - ✅ **C12** `scanner_parity.rs` — **copied, as the row said; the 1 test kept.** Production's
    post-filter (`ssr::catalog::postprocess::build_loci` and its helpers, 220 lines) is now a
    `production_post_filter` module inside the test, unchanged but for taking the golden
    catalog's four settings, writing `Motif::new`'s length check inline, and returning a plain
    locus. Before production's copy was removed, the two were run side by side on the scanner's
    intervals over the synthetic reference and returned the same 17 loci. The test's readout is
    unchanged: 16 of 16 golden loci recovered, 15 exact, 1 boundary wobble, 1 scanner-only.
    **Deviation, recorded 2026-09-13, which changes C19's action:** the golden catalog is a
    16-locus BGZF text file, so rather than freeze its contents, ng gained a 50-line test-only
    reader for it, `ng::golden_catalog`, which C19's three tests will use too. The committed file
    stays the oracle.
  - ✅ **C13** `window_coverage/production_parity.rs` — **frozen as digests; the 1 test kept.**
    Production's accumulator emitted 447,581 windows over the 200 generated streams, about 13 MB
    window by window, over the fixture cap. What it emitted on each stream was recorded at
    `d9e7b076` as a count and an FNV-1a digest of every window's contig, centre, mean depth and GC
    fraction (`window_coverage/testdata/production_windows_by_seed.tsv`, 200 rows), and ng's stream
    is reduced the same way. **Given up:** a divergence still fails at its seed, but no longer names
    the window or the field. Changing one digest fails the test.
- ✅ **C14. The four parity examples** — `ng_psp_against_production.rs`, `ng_psp_parity.rs`,
  `ng_psp_head_encoding.rs` (`test = true` in `Cargo.toml`), `paralog_score_parity.rs` — are
  deleted with a line each in the report saying which document already holds their result.
  *Depends:* — · *Source:* PROJECT_STATUS entries for each.
  Done 2026-09-13. Where each result lives:
  - `ng_psp_against_production.rs` — the reader's speed against production's `.psp` (1.8× on 62
    samples): `psp_file_format.md` §5.4 and `reports/reviews/perf_ng-psp_2026-08-30.md`.
  - `ng_psp_parity.rs` — ng's psp round-trips production's records and agrees with an independent
    encoder: `reports/implementations/ng_psp_h1_2026-08-28.md` and `…_h4_2026-08-30.md`.
  - `ng_psp_head_encoding.rs` — variable-length against fixed-width record heads after compression:
    `psp_file_format.md` §4.3 and `reports/implementations/ng_psp_head_h2_2026-09-04.md`,
    `…_h3_2026-09-04.md`.
  - `paralog_score_parity.rs` — the hidden-duplication filter's statistics on tomato2 data:
    `reports/implementations/paralog_r1_data_validation_2026-07-01.md`. Its Python companion
    `benchmarks/tomato2/src/paralog_score_parity.py` stays with the benchmark's other scripts
    (plan §8, benchmark drivers are D6's).
  The two with `test = true` lose their `[[example]]` entries in `Cargo.toml` — and with them the
  36 tests those two harnesses carried on their own mechanics (4 and 32), which is why the suite
  drops by 36 here — and four doc
  references that pointed at them as live tools now say they were deleted.

#### The rest of C, found by sweeping for oracles outside the thirteen dedicated files

**Recorded 2026-09-12, before C1 was started.** §1 said "thirteen files, 103 tests"; that counted
only the files whose *name* says parity. Sweeping `src/ng/` and `src/pop_var_caller_exp/` for every
site that names a production module found **25 further files**. Measured on the tree at `d9e7b076`:

```
grep -rnE 'crate::(pileup|psp|pileup_record|genetics|pop_var_caller|var_calling|vcf|ssr|
           paralog|sample_summary|baq|norm_seqs)(::|;| )' src/ng src/pop_var_caller_exp src/main_exp.rs
```

128 lines, of which 13 are doc comments. The thirteen dedicated files hold 103 tests over 11,927
lines; the twenty-five further files hold **46 more tests that execute production code**, spread
over 18 of them. The other seven are not oracles at all, and they are the reason this sweep
mattered:

- **Three of them are shipped code, not tests.** `src/pop_var_caller_exp/` — ng's own command
  surface, which §1 counts as part of ng — still imports `crate::pop_var_caller::common`. B9's
  guard is scoped to `src/ng` and cannot see them, so Checkpoint B passed with them present. They
  must land before D1 deletes `src/pop_var_caller/`.
- **One is a test-fixture builder that shared infrastructure also uses.** `cram_files` writes
  synthetic FASTA and CRAM files to a tempdir. It lives in `src/pileup/per_sample/`, which D3
  deletes, and `src/bam/` (three files) and `src/fasta/` (one) call it from their own tests. Those
  are modules §2 keeps, so it is rehomed, not frozen.
- **Two are B3 repoints that were listed and not made.** B3's row named
  `benches/ng_psp_perf.rs` and `examples/dhat_ng_psp.rs`; both still say
  `pop_var_caller::pileup_record::ChainId`.

| step | file(s) | tests | what it runs from production | action |
|---|---|---:|---|---|
| C15 | `pop_var_caller_exp/{generate_psps,typed_regions,calling_run}.rs` | — | `pop_var_caller::common::{current_command_line, rfc3339_now, DEFAULT_BUFFERED_IO_CAPACITY}` | **copy into ng.** Shipped code — **own commit**, and it is Milestone B's work arriving late |
| C16 | `pileup/per_sample/cram_files.rs` and its nine callers | — | the synthetic FASTA/CRAM builder itself | **rehome to `src/bam/`** — see the deviation note below |
| C17 | `benches/ng_psp_perf.rs`, `examples/dhat_ng_psp.rs` | — | `pileup_record::ChainId` | **repoint** to `ng::types::ChainId`, finishing B3 |
| C18 | `alignment/stutter.rs`, `calling/likelihood/mod.rs`, `locus_generation/pileup/{generator,mod}.rs`, `psp/{mod,header}.rs` | 6 | six constants asserted equal to ng's copies: `MAX_SLIP`, `MIN_BASE_ERROR`, five walker `DEFAULT_*`, `psp::header::HEAD_MAGIC` | **freeze** — the fixture is the literal value, with the production commit that wrote it |
| C19 | `repeat_catalog/anchor.rs`, `region_typing/mod.rs`, `reference_info.rs` | 3 | `ssr::catalog::io::CatalogReader`, reading the committed `tests/data/tandem_repeat/golden.ssr_catalog.bed.gz` | ~~freeze the file's contents~~ **repoint to `ng::golden_catalog`**, the reader C12 added — the golden catalog stays the oracle without production's parser |
| C20 | `calling/likelihood/generic.rs` | 1 | `var_calling::per_group_merger::standard_log_likelihood`, `pileup_record::AlleleSupportStats` | **freeze** |
| C21 | `calling/genotype_prior/dirichlet_multinomial.rs` | 1 | `genetics::dirichlet_multinomial_log_priors` | **freeze** |
| C22 | `alignment/ssr_marginal_sequence.rs` | 3 | `ssr::cohort::pair_hmm::{HmmScratch, align_subst}` | **freeze** |
| C23 | `locus_generation/ssr.rs` | 6 | `ssr::pileup::{fetch_reads, alignment, footprint, locus_tally}`, `ssr::types` | **freeze** |
| C24 | `region_typing/segment_criteria.rs` | 9 | `ssr::catalog::postprocess::build_loci`, `ssr::catalog::{CatalogParams, trf::TrfRecord}`, `ssr::types::{Locus, Motif}` | **freeze.** A further 21 tests in the file only reach production through the helper that builds *both* sides' settings from one source; those need the helper split, not a fixture |
| C25 | `ref_seq.rs` | 1 | `pileup::per_sample::read_processor::RawContigRefCache` | **freeze** |
| C26 | `read/prepared_read.rs` | 6 | `pileup::walker::{PreparedRead, MateRole}` through the four `#[cfg(test)]` bridges B7 left in place | **freeze**, and delete the bridges. *Depends:* C1 and C9, which consume them |

- ✅ **C15** · **own commit** — done 2026-09-13. `current_command_line`, `rfc3339_now` and the
  civil-date helper it calls, with their two tests, moved unchanged into a new
  `pop_var_caller_exp::provenance`; `DEFAULT_BUFFERED_IO_CAPACITY` (64 KiB), used once, became a
  constant in `typed_regions.rs`. Three doc links in the command surface that pointed at production's
  CLI module became plain text, since they would break `cargo doc` once it is deleted. The identity
  oracle cannot see this move — it strips the `##commandline` line and the psp header's timestamp,
  the two things these functions write — so what stands behind it is that the code is copied
  byte for byte and its two tests pass. · *Depends:* B9 · *Source:* §1 (ng includes `pop_var_caller_exp`).
- ☐ **C16** · *Depends:* — · *Source:* the B5 deviation, which settled the same question for `CigarOp`.
- ☐ **C17** · *Depends:* — · *Source:* B3's row.
- ☐ **C18**–**C25**, one checkbox each, in table order · *Depends:* B9.
- ☐ **C26** · *Depends:* C1, C9.

  **Deviation, recorded 2026-09-12: `cram_files` goes to `src/bam/`, and §2's "shared
  infrastructure is not touched" gives way.** The alternative is a copy in ng and a second copy
  left behind for `bam` and `fasta` to use — two builders of the same synthetic CRAM, drifting.
  §2 excluded shared infrastructure to keep this plan from *pruning* it; taking in a file whose
  four non-ng callers already live there is not pruning, and B5 already crossed this line for
  `CigarOp` for the same reason. Production's own `ssr` callers keep compiling through the move
  because the new path is a `use`, and they go at D4 regardless.

> **Checkpoint C: production is now imported only by production.** The sweep command at the head
> of the second table, run over `src/ng`, `src/pop_var_caller_exp`, `src/main_exp.rs`, `tests/ng_*`,
> `benches/ng_*` and `examples/ng_*`, returns nothing; `cargo test` runs the frozen tests green;
> A2's md5s equal A1's.
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
| C — freeze | md5s = A1; frozen tests green; the sweep at the head of §5.C's second table, run over ng's source, command surface and `ng_*` tests/benches/examples, returns nothing |
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
