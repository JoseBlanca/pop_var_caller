# Code Review: read_phase
**Date:** 2026-11-03
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the read-backed phase build, branch `neighbour-pass`, diff `eedd7a63..a5941dc1`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `eedd7a63..a5941dc1`, eleven commits that build
  [`neighbour_evidence_read_phase.md`](../../ng/spec/neighbour_evidence_read_phase.md): the
  `refine-genotypes` removal, psp format 1.1's chain-id base, the capture in the merge, read alleles
  through calling, the fragment keeper, fragments through the writer and the spill, the link
  log-odds, phase from reads alone, reads in the neighbour-evidence pass, and the invariance tests.
- **Reviewed against:** commit `a5941dc1`.
- **In-scope files:** `src/read_phase/*`; `src/neighbour_evidence/{phase,refine,stage}.rs`;
  `src/run/cohort_merge/{build,rounds,serial}.rs`; `src/run/{records,callers,psp_caller,edge_deletion_gathering,report,psp_source}.rs`;
  `src/run/paralog_filter/{pass_one,pass_three,spill,spill_file}.rs`; `src/vcf/{writer,header,assemble,mod}.rs`;
  `src/psp/{record,header,mod}.rs`; `src/locus_generation/mod.rs` and `pileup/{chain_id_allocator,generator,genome_walk,active_read_set}.rs`;
  `src/cli/{calling_run,call_from_psps,call_from_alignments,command_line,mod,test_fixtures,mode_equivalence}.rs`;
  `src/main.rs`; `examples/ng_str_stutter_benchmark.rs`; the tests beside each.
- **Deliberately out of scope:** files where the diff only adds a zero `chain_id_base` or an empty
  `read_alleles`/`fragments` to a test fixture; the population phase output (reviewed 2026-11-02);
  `vendor/`; documentation prose (the report and user pages came after `a5941dc1`).
- **Categories dispatched**, in six reviewers of two categories each, every one in its own git
  worktree checked out at `a5941dc1`: reliability and errors (always); naming and defaults (new
  constants, an option, header lines); idiomatic and smells (always); refactor_safety and
  module_structure (multi-file, new module); unsafe_concurrency and float_portability (rayon, a
  `OnceLock`, log-odds sums that decide output); tooling and extras (byte decoders, stable output,
  hot paths, commit messages against the diff). Six rather than twelve to stay within the session's
  agent budget; each reviewer applied both checklists in full. Per-category files:
  `tmp/review_2026-11-03_read-phase/`.

## 2. Verdict

**Approve-with-changes.** No blocker. Eight major findings, none of which changes what the build
writes on any input the tests or the measurements produced; they are a possible panic on a record
pair the edge-deletion gathering is expected never to leave, two places where the type system lets
phase be dropped silently, an option whose effect is not stated where it is set, and four gaps in
the tests — chiefly that nothing end to end would notice the read cap losing the walked region's
first chain id. The arithmetic was checked portable and thread-invariant, including a temporary
test of the pass's read-backed link test at 1 and 8 threads.

## 3. Execution status

Run by the orchestrator at `a5941dc1`, release profile, native cargo on macOS arm64:

- `cargo fmt --check` — exit 0, no diff.
- `cargo clippy --release --all-targets --all-features -- -D warnings` — exit 101, from
  `tests/ng_calling_loop_allocation.rs` only (built only with the `dhat-heap` feature; `missing field
  low_quality_reads in initializer of pop_var_caller::calling::GenericObservation`; the file was last
  changed by `bbc757c9`, before this branch). Without `--all-features`: clean.
- `cargo test --release --lib` — 5,121 passed, 2 failed (`an_offset_past_the_reference_panics`,
  `lgamma_panics_on_non_positive_in_debug`: `should_panic` tests of debug-only assertions, which a
  release build never fires; pre-existing). Integration tests: 33 passed.
- `cargo doc --release --no-deps` — exit 0, no warnings from this crate.
- Not run: `cargo audit` (not installed).
- Reviewers' own runs: the reliability reviewer ran 14 mutations (2 killed, 10 survived with a changed
  behaviour shown by a probe test, 1 changed no reachable behaviour, 1 survived with no change shown);
  the tooling reviewer re-ran the `chain_id_base → 0` mutation against the whole suite.
- Findings labelled "Needs verification": 1 (M1, Medium confidence: the panic's input was not built).

## 4. Open questions and assumptions

1. **Should phase from reads alone have its own switch?** (M7.) The spec rules that every diploid run
   without refinement writes it; it does not say whether a user can turn it off. Recommended: no new
   switch for now, but say in `--refine-genotypes`'s help what happens where refinement does not run.
2. **Can the edge-deletion gathering leave a left-padded generic record and an STR record at one
   written position unmerged?** (M1.) Its own counter for that case is expected to stay zero; if it
   can, the keeper panics.

## 5. Top 3 priorities

1. **M5 — a test that the walked region's first id reaches the cap.** Two reviewers independently
   replaced it with zero and the whole suite stayed green; it is the property psp format 1.1 exists for.
2. **M2 — fragments off `SpillEntry`, beside it by type.** One wrong reader call in pass three would
   drop all read-backed phase from every run with the duplication filter on, silently.
3. **M1 — the keeper's order check.** A panic, however unlikely, on a stream the writer accepts;
   the fix is one line.

## 6. Findings

### Blocker

None.

### Major

**M1: [src/read_phase/keeper.rs](../../../../src/read_phase/keeper.rs#L197) — the order check uses
`region.start`; the stream is ordered by written `POS`.** *Confidence: Medium (needs verification).*
The writer admits one tie: a generic record then an STR record at the same written `POS`. If the
generic record is padded left (its span starts at P+1, written at P) and the STR record starts at P,
the keeper's `assert!(newest.at <= at)` fails and the run panics. That needs the edge-deletion
gathering to leave the pair unmerged, which its `generic_records_left_unmerged` counter expects never
to happen. **Fix:** key the keeper on the written position (`place_of`'s), or replace the assertion
with ordering on `max(previous, start)`; add a test of that tie. *Category: reliability.*

**M2: [src/run/paralog_filter/spill_file.rs](../../../../src/run/paralog_filter/spill_file.rs#L255) —
pass two's reader and pass three's return the same `SpillEntry`.** *Confidence: High.* Its
`fragments` is empty from `read()` and filled from `read_with_fragments()`; pass three switched to
`read()` would compile and write every line without fragments, the trap spec §6 names. **Fix:** take
`fragments` off `SpillEntry`; encode them as a separate argument; give pass three a reader whose item
is `(SpillEntry, RecordFragments)`. *Category: refactor_safety.*

**M3: [src/read_phase/stage/tests.rs](../../../../src/read_phase/stage/tests.rs#L39) — no stage test
runs a line with fragments through the streaming path.** *Confidence: High.* Every fixture spans at
most 400 bases, so every line settles in `finish()`; forgetting fragments after 1,000 bases instead of
2,000 passes the suite (mutation shown). **Fix:** a fixture of linked records spread over more than
3,000 bases, asserting the same output as when pushed whole. *Category: reliability.*

**M4: [src/read_phase/stage/tests.rs](../../../../src/read_phase/stage/tests.rs#L134) — "skipped, not
a gap" cannot tell the phaseability rule from "any fragment".** *Confidence: High.* Its middle record
has no fragments; a middle record whose one fragment is seen nowhere else would also pass if the rule
were weakened. **Fix:** give the middle record a fragment at no other record. *Category: reliability.*

**M5: [src/cli/mode_equivalence.rs](../../../../src/cli/mode_equivalence.rs#L560) and
[src/run/cohort_merge/build.rs](../../../../src/run/cohort_merge/build.rs#L771) — nothing end to end
checks that the walked region's first id reaches the cap, or that read alleles are renumbered through
`AlleleRemap`.** *Confidence: High.* *Categories: reliability, tooling/extras (convergent).* Replacing
the base with 0 at `build.rs:771`, and skipping the remap in `read_alleles_on_the_record`, each
passed all 5,131 lib tests. The deep invariance test cannot fail: every read carries the pair in cis,
so any 32 kept give the same phase. Spec §7.5's direct comparison of captured `(rank, allele)` lists
is not built. **Fix:** compare `LocusReads` captures of the second walked region between a whole-contig
and a `--regions` walk as `(id − base, allele)` lists of length 32; give the deep fixture a minority of
trans reads so which 32 are kept changes the output; a fixture where selection drops an allele.

**M6: `src/read_phase/{record_reads.rs:48, keeper.rs:210, keeper.rs:249, keeper.rs:286, stage.rs:197,
stage.rs:281, stage.rs:353}` — seven `expect` calls or map indexings without a `// PANIC-FREE:`
comment.** *Confidence: High.* The project's rule. **Fix:** a comment naming the invariant at each, or
return an error where the invariant is the caller's. *Category: errors.*

**M7: [src/vcf/writer.rs](../../../../src/vcf/writer.rs#L153) — phase from reads alone is on in every
diploid run without refinement, and the command line does not say so.** *Confidence: Medium.* Ruled
by the spec; the gap is that `--refine-genotypes`'s help describes refinement only, so a user who sets
`off` gets `|` and `PS` without having asked. (The user pages committed after `a5941dc1` now describe
it.) **Fix:** one sentence in the option's help on both commands; open question 1. *Category: defaults.*

**M8: [src/read_phase/stage.rs](../../../../src/read_phase/stage.rs#L437) — `with_phase` assumes `GT`
is the first key and appends `PS` unconditionally.** *Confidence: Medium.* `eligible_columns` looks
`GT` up by position while the rewrite writes the first field; and a `FORMAT` already ending in `PS`
would get two. Neither arises from calling's own lines today. **Fix:** rewrite at `GT`'s position and
reuse an existing `PS` key, as `neighbour_evidence::line::rewrite_line` does. *Categories:
idiomatic/smells (cross-category), refactor_safety.*

### Minor

Each in full, with its fix, in the per-category file named.

- Mi1 `src/run/records.rs:603` — `zip` pairs read-allele lists with samples and truncates a mismatch silently. *reliability, refactor_safety (convergent).*
- Mi2 `src/run/paralog_filter/spill.rs:606` — the fragment decoder has no malformed-input test, and leftover bytes are reported as an out-of-range length. *tooling, reliability.*
- Mi3 `src/psp/record.rs:2509` — a damaged `chain-id-base-change` decodes silently and carries to every later record of the block. *tooling.*
- Mi4 `src/run/paralog_filter/spill/tests.rs:1053` — the spill's property test never generates fragments. *reliability.*
- Mi5 `src/neighbour_evidence/stage/tests.rs:104` — the pass's thread-count tests pass no fragments, so the read-backed link test never runs in parallel (a temporary test was identical at 1 and 8 threads). *concurrency.*
- Mi6 link-test fixtures use a copying length of 300,000 bases, so five guards in `phase_heterozygotes`/`attach_by_reads` can be removed without a failure. *reliability* (filed Major there; Minor here because each guard is also a spec rule a reader can check, and none changes output on the measured cohorts).
- Mi7 `src/read_phase/link.rs` / `src/neighbour_evidence/phase.rs` — allele indices are `u32`, `u16` or `usize` by module, and `SiteReads::link` narrows with `as u16`. *naming, idiomatic, reliability (convergent).*
- Mi8 `src/vcf/header.rs:763` — `##readBackedPhase=` omits the break threshold the stage applies; the cap of 32 is in no header line. *defaults.*
- Mi9 `src/cli/calling_run.rs:509` — the read-check log hard-codes 199 to 1, GQ 30 and 1 kb. *defaults, refactor_safety (convergent).*
- Mi10 `src/run/records.rs:773` — `evidence_for_output` fills `read_alleles` with a default both drivers overwrite; a driver that forgets writes no read alleles. *refactor_safety.*
- Mi11 `src/neighbour_evidence/stage.rs:200` — fragments kept beside sites in hand-synchronised parallel vectors. *refactor_safety.*
- Mi12 `src/neighbour_evidence/phase.rs:560` — `AddAssign for PassReadCounts` without a destructure. *refactor_safety.*
- Mi13 `src/locus_generation/pileup/open_record.rs:1227` — minted with base 0 and corrected in the generator; a path that bypasses it keys the cap on raw ids. *refactor_safety.*
- Mi14 `src/psp/record.rs:610` — the legacy-layout ladder is unreachable since format 1.1, and its oldest rung would read a zero base. *refactor_safety.*
- Mi15 `src/read_phase/keeper.rs:17`, `capture.rs:5` — `read_phase` and `vcf` import each other, and `read_phase` reaches into `run::psp_source::read_cap` for its hash. *module_structure.*
- Mi16 `src/read_phase/stage.rs:366` — the reads-only stage parses and allocates per line and per column on every diploid run; its parsing duplicates the pass's line reader. *idiomatic/smells.*
- Mi17 `src/vcf/writer.rs:95` — two `Option` stage fields that must never both be `Some`. *idiomatic/smells.*
- Mi18 `src/neighbour_evidence/phase.rs:349` — `phase_heterozygotes` is 137 lines; the relation/orientation arithmetic is written out three times; attached records carry a fabricated `Candidate`. *idiomatic/smells.*
- Mi19 `src/run/cohort_merge/build.rs:773`, `spill.rs:530`, `keeper.rs:180` — per-sample clones and per-record vectors on paths every record takes. *idiomatic/smells.*
- Mi20 naming: positional tuples in the keeper, adjective type names in the stage, `reads_phase_alone`, `qualities` for genotype qualities, bare VCF column numbers. *naming.*
- Mi21 `src/read_phase/keeper.rs:93` — the keeper's and the stage's cost is measured only at 50 and 200 samples, not at the 3,000 the spec's §4.2 worries about. *tooling* (the implementation report gives 0.04% for the keeper at 200).
- Mi22 the stage writes an unreadable line without `PS` and logs nothing; no test that a contig change clears the stage's counts and open block. *reliability.*

### Nits

Grouped in the per-category files (23 in all): doc-comment wording, a broken `write_line` doc
reference, test names, constant placement.

## 7. Out of scope observations

- `tests/ng_calling_loop_allocation.rs` does not compile under `--all-features` (pre-existing, from
  `bbc757c9`).
- Two `should_panic` tests of debug-only assertions fail in every release-build test run.

## 8. Missing tests to add now

From the reliability file's `## Missing tests`, with code there:

- `the_cap_keeps_the_same_ranks_whole_contig_and_regions` (M5) — `LocusReads` captures compared as `(id − base, allele)` across a whole-contig and a `--regions` walk.
- `the_stage_phases_alike_when_lines_stream_past_two_reaches` (M3).
- `a_record_whose_only_fragment_is_seen_nowhere_else_is_skipped` (M4).
- `a_read_made_link_beyond_the_copying_length_is_not_made` and the four sibling guards (Mi6).
- `a_padded_generic_record_and_an_str_record_at_one_position_pass_the_keeper` (M1).
- `a_damaged_fragments_section_is_refused` and `a_damaged_chain_id_base_change_is_refused` (Mi2, Mi3).
- `the_pass_phases_alike_with_reads_at_one_thread_and_eight` (Mi5).

## 9. What's good

- The psp head carries the base as a signed change from the previous record, restarting per block, so
  a reader that skips bodies still rebuilds it — and a round-trip test under every skip pattern pins it
  ([src/psp/record.rs](../../../../src/psp/record.rs)).
- The weight table is computed once from `crate::float`, the sum runs in fragment order, and the
  ceiling at one half makes a Phred-0 read weigh exactly zero, tested on hand values
  ([src/read_phase/link.rs](../../../../src/read_phase/link.rs)).
- The capture reuses the merge's own allele derivation (`assemble_capturing`), so a read's allele is
  the one its support was tallied under, rather than a second derivation that could drift
  ([src/run/cohort_merge/build.rs](../../../../src/run/cohort_merge/build.rs)).
- The byte digests were re-recorded only after showing that removing the new lines reproduces the old
  checksums exactly ([src/cli/cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs)).

## 10. Commands to re-verify

- `cargo fmt --check`; `cargo clippy --release --all-targets -- -D warnings`;
  `cargo test --release --lib`; `cargo test --release --no-fail-fast --test '*'`.
- `cargo test --release --lib mode_equivalence read_phase neighbour_evidence::phase paralog_filter::spill`.
- The M5 mutation: replace `records.first().map_or(0, |record| record.chain_id_base)` with `0` at
  `src/run/cohort_merge/build.rs:771`; the new test from section 8 must fail.
