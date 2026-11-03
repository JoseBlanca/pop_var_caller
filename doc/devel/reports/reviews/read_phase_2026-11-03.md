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

## 11. Dispositions

Applied 2026-11-03 on branch `neighbour-pass`, commits `dba90cca`..`1e531539`. Every new test was
run against the mutation it is meant to catch and failed; the mutations are named in the commit
messages. The full lib suite passes apart from the two pre-existing debug-only `should_panic`
tests; clippy is clean apart from the two vendored noodles-cram warnings.

| Finding | Disposition |
|---|---|
| M1 keeper ordered by `region.start` | Fixed in `54a8fef6`: keyed on `written_position`, the writer's own key. The tie test panicked on the old key. |
| M2 one `SpillEntry` for both readers | Fixed in `a8ea1016`: fragments are not a field; pass three's reader is `SpillReaderWithFragments`, item `(SpillEntry, RecordFragments)`. |
| M3 no streaming stage test | Fixed in `30934ba3`: two tests past `2R`; the first kills forgetting after `R` (s3). |
| M4 "skipped, not a gap" too weak | Fixed in `30934ba3`: the middle record has a fragment seen nowhere else; kills s1. |
| M5 base and remap untested end to end | Fixed in `40a47878`: whole-contig vs `--regions` psps, the merge's capture compared as `(id − base, allele)`, 32 entries; plus two unit tests. Kills c1, c2 and a dropped generator stamp. The minority-of-trans-reads fixture was not built: one fragment weighs about 3.9 nats, so a minority of trans reads cannot change a 32-read link and the VCF test still could not fail. |
| M6 `expect` without `PANIC-FREE` | Fixed in `54a8fef6`, `30934ba3`; the keeper's map index became `get` with an `UNREACHABLE` comment naming the forgetting rule. |
| M7 option help silent on reads-alone phase | Fixed in `1e531539`: help on both commands and the user page. No new switch. |
| M8 `with_phase` writes the first field, appends `PS` twice | Fixed in `30934ba3`: `GT` at its `FORMAT` place, an existing `PS` reused; test kills "GT first". |
| Mi1 `zip` pairing | Fixed in `40a47878`: the capture names each entry's sample. |
| Mi2 fragment decoder unchecked | Fixed in `a8ea1016`: sample and Phred range checks, `SectionLongerThanItsEntries`; one test per refusal, cuts with both readers. |
| Mi3 damaged `chain-id-base-change` | Fixed in `16a441e9`: a base above the record's smallest read is refused. Covers generic records; STR and bundle records carry base 0 and are not covered; damage that leaves the base at or below the smallest read still passes. |
| Mi4 property test without fragments | Fixed in `a8ea1016`. |
| Mi5 thread test without reads | Fixed in `a07519e0`: 1 and 8 threads with mixed fragments, identical, and different from no reads. |
| Mi6 five guards untested | Fixed in `a07519e0`: one test each, kills p1, p2, p3, p4, p6; and r1 (link-test `GQ` from the refined posterior), which survived the whole suite. |
| Mi7 allele index types | Partly: every narrowing is now checked (`try_from`) and `Heterozygote` can only be built by `of`. **Deferred:** an allele-index newtype. `AlleleId(u16)` already exists and is the right type, but adopting it in four structs and every fixture changes no behaviour now that the narrowing is checked. |
| Mi8 header lacks `break_min`, cap invisible | Fixed in `1e531539`: both header lines, spec §3, and a logged count of what the cap removed. The three calls digests re-recorded; with the old header value they reproduce exactly. |
| Mi9 log hard-codes thresholds | Fixed in `1e531539`. |
| Mi10 placeholder `read_alleles` field | Fixed in `40a47878`: an argument of `assemble_record`. |
| Mi11 parallel vectors; empty = no reads | Fixed in `a07519e0`: `HeldSites`; `Option<&[RecordFragments]>`. |
| Mi12 `AddAssign` without destructure | Fixed in `a07519e0`. |
| Mi13 base minted as 0, stamped later | Kept as built: every generic record leaves through the one stamping branch, and the end-to-end test of M5 now fails if the stamp is dropped (shown). |
| Mi14 unreachable legacy ladder | Fixed in `8b597ed1`: about 320 lines removed; a manifest without `chain-id-base-change` is refused, tested. |
| Mi15 `read_phase` ↔ `vcf`/`run` cycle | Fixed in `dba90cca`: `read_phase::fragments`, `run::fragment_keeper`, `crate::splitmix`. `read_sampling`'s finaliser is a different function (no golden-gamma add) and was left alone. |
| Mi16 per-line allocation, duplicated parsing | Fixed in `30934ba3`: `vcf::line_fields` shared by the stage and the pass's reader; no per-column `Vec`. |
| Mi17 two optional stages | Fixed in `1e531539`: `LineStage`. |
| Mi18 long functions, fabricated `Candidate` | Fixed in `30934ba3` (stage `Join`) and `a07519e0` (`link_to_next`, `written_phases`, `Member` without `Candidate`, `Relation::between`/`carry`). |
| Mi19 per-record allocations | Fixed: keeper buffer (`54a8fef6`), capture (`40a47878`), spill (`a8ea1016`), stage and `write_line` by value (`30934ba3`, `1e531539`); `attach_by_reads` trees replaced by vectors (`a07519e0`). |
| Mi20 naming | Fixed across the groups (tuples to structs, `BlockPlacement`, `BlockRecord`, `SampleReadBlocks`, `upper_on_copy_0`, `genotype_qualities`, `writes_phase_from_reads_alone`, `read_phase_tally`). |
| Mi21 cost at 3,000 samples unmeasured | **Deferred:** needs a 3,000-sample simulated cohort; the implementation report's 200-sample figure (keeper 0.04% of the run) is the measurement spec §4.2 asks for first. None of these fixes adds per-record work beyond one counter. |
| Mi22 unread line silent; contig change untested | Fixed in `30934ba3`: counted and warned; contig test kills dropping the per-contig clear. |
| `write_stream` writes no fragments | Fixed in `1e531539`: items are `(record, fragments)`. |
| Commit `31669bfd` numbers do not reproduce | Cannot amend history. The reviewer's run at `a5941dc1` is the figure to use: 39 records written, 333 read alleles reaching the keeper, 222 kept. |
| Nits: doc splits in `# Errors`, wrong capture-test comment, `pub` docs, `cap_key` visibility, `MIN_FRAGMENT_ENTRY_BYTES`, Phred/floor coupling test, `num_obs` rule, `error_phred` only when capturing, bool parameters, struct-update reader | Fixed in the groups above. |
| Nit: `chain_id_base` typed `ChainId` | Skipped: `ChainId` is an alias of `u64`, so the change gains no checking. |
| Nit: `Arc<str>` clone per paired read in the allocator | Skipped: the reviewer judged it acceptable beside BAM decoding; unmeasured. |
| Nit: three `refine_window_*` wrappers; `OrderedRecordStages` | Deferred: the reviewers mark both optional, and the wrappers now pass an explicit `None`. |
| Nit: `KeeperCounts` `pub(crate)` | Skipped: it is the return type of the keeper's public `counts()`. |
| Cross-category: link-test `GQ` taken before renormalising | No change: renormalising can only raise `GQ`, so this can refuse a link, never make a wrong one; now tested (r1). |
