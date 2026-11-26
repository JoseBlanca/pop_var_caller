# Code-quality map: input/output layer — `psp`, `read`, `bam`, `alignment`, `read_phase` (2026-11-26)

Read-only review. Nothing was built or run. About 70,000 lines are in scope; `src/read`, `src/bam`
and `src/alignment` + `src/read_phase` were each read by a delegated reviewer, and the findings
below that carry call-changing consequences were re-checked against the code before being
written here.

## Summary, in priority order

1. **A BAM cut off at a compressed-block boundary reads as a shorter, valid file.** Every region
   after the cut returns zero reads and the run finishes without an error (bug B1, likely).
2. **A record with bases but no base qualities (`QUAL *`) stops the whole run** instead of being
   dropped and counted (bug B3, likely).
3. **Index lookup misses `sample.bai` and never notices a stale index.** A Picard/GATK-style
   index beside the file is reported missing, and an old `.csi` is preferred over a fresh `.bai`
   (bugs B4, B5).
4. **CRAM: a long read in an earlier container can be skipped when a cursor jumps to a region**,
   so CRAM and BAM of the same data can disagree (bug B2, suspected; long or spliced reads only).
5. **About 3,500 of `src/bam`'s 5,600 lines are dead**: an older reader stack hidden by
   `allow(dead_code)`, duplicating the live one in `src/read/input`. Delete it, then fold what
   remains of `src/bam` into `src/read` (A1, A2).
6. **`read` → `psp` dependency:** all four uses are one function that turns the read-filter
   settings into psp-header provenance. Move it to `run`, beside the code that reads those keys
   back; then neither `read` nor `psp` depends on the other (A3).
7. **The psp encoder and decoder do share one description of the format — as a checked
   fingerprint, not as a driver.** Each section's encoder and decoder sit in the same file, the
   field table is written into every file's header and checked on open, and a golden-bytes test
   keeps the hand-written codec in step. That is a sound design; leave it (A5). What needs work is
   psp's documentation, which still describes deleted code and renumbered format versions (R1).
8. **A soft-clipped read keeps an edge indel that the same read unclipped would drop** (B6,
   suspected; changes output, owner decision).

---

## Bugs

### B1 — A BAM truncated at a BGZF block boundary is read as complete, losing every later region

- **Where:** `src/read/input/aligned_reads_reader/bgzf_blocks.rs:270-275` (`read_frame_into`),
  `:251-252` (`inflate_at`), `src/read/input/aligned_reads_reader/bam.rs:233-241`.
- **What:** `read_frame_into` treats `UnexpectedEof` while reading the 18-byte block header as a
  clean end of stream, and `inflate_at` then returns an empty block. noodles reports
  `bytes_read == 0`, and `bam.rs` reads that as "end of file inside a chunk — try the next one".
  When the index was built from the complete file, every chunk past the cut seeks beyond the end,
  gets an empty block, and is skipped. **Scenario:** a BAM copied with an interrupted `cp`/`rsync`,
  or written by a killed writer, whose last byte falls between blocks (or inside the first 18
  bytes of one). Every read after that point, to the end of the file, is missing; the sample's
  calls there become no-calls or reference calls, and nothing is logged. No check for the 28-byte
  BGZF end-of-file marker exists (grep finds none in `src/read` or `src/bam`). The existing test
  `t10_a_truncated_file_fails_once…` (`src/read/input/open_bam.rs:2005`) cuts at three quarters of
  the file size, which lands mid-block and so exercises only the case that already errors.
- **Proposal:** at open, check the file ends with the BGZF EOF marker; refuse (or at least warn)
  when it does not. Treat a block header cut short as `InvalidData`, not end of stream. Add a
  truncate-at-block-boundary test. The CRAM side makes the same assumption
  (`aligned_reads_reader/container.rs:552-570`: zero bytes from `read_container` is end of
  stream); checking for the CRAM EOF container at open covers it.
- **Downside:** a file genuinely lacking the marker (very old writers) would be refused; a warning
  is the softer choice.
- **Effect:** output changes only on damaged input (a run that silently lost data now stops).
  Speed/memory none. **Cost:** S. **Confidence:** high (likely, traced).

### B2 — CRAM positioning can skip an earlier container whose long read reaches the region

- **Where:** `src/read/input/aligned_reads_reader/cram.rs:344-363` (`first_entry_reaching`).
- **What:** the function steps back over `.crai` entries and stops at the first one whose span does
  not reach `start`. An entry before that can still reach it. **Scenario:** slice S1 starts at
  1,000 with span 50,000 (one read with a 49 kb `N` skip, or a long read with a big deletion);
  slice S2 covers 1,100–1,400. A jump to 20,000 stops at S2, and the S1 read is not seen for that
  region. It fires only on a jump (a cursor's first region, or a backward move), so whether a read
  is lost can depend on where worker boundaries fall. The BAM path uses index bins and has no such
  gap.
- **Proposal:** when the `.crai` is grouped by contig (`src/read/input/open_bam.rs:727`), store
  the running maximum of `start + span`; `first_entry_reaching` becomes one `partition_point` on it.
- **Downside:** 8 bytes per slice; possibly a few more container decodes after a jump.
- **Effect:** output possible (rare inputs). **Cost:** S. **Confidence:** medium (suspected; no
  fixture).

### B3 — A read with bases but no base qualities stops the run

- **Where:** `src/read/aligned_read.rs:159` (copies an empty `quality_scores()` as is);
  `src/read/prepared_read.rs:350-356` (`length()` refuses `seq.len() != bq_baq.len()`);
  `src/locus_generation/pileup/genome_walk.rs:1024` (propagates that as a `WalkerError`).
- **What:** SAM allows `QUAL *` with a real `SEQ`. The mismatch filter tolerates it
  (`src/bam/alignment_input.rs:723-732`); nothing drops the read, so it reaches the generic-path
  walker, where `length()` fails and the walk returns the error. **Scenario:** a BAM produced from
  FASTA input, or some assembly/long-read BAMs — one such read aborts the sample. A record with
  `SEQ *` reaches the same failure when `min_read_length` is off.
- **Proposal:** a counted drop reason in the read filter ("no base qualities"), or a documented
  default quality filled at decode. The first is simpler and visible in the tallies.
- **Downside:** one new tally field and one provenance key.
- **Effect:** output yes (a run that aborted now finishes). **Cost:** S. **Confidence:**
  medium-high — I did not confirm every caller of the walker treats this error as fatal.

### B4 — `sample.bai` / `sample.crai` (extension replaced, not appended) are not found

- **Where:** `src/bam/index_preflight.rs:306-317` (`existing_index_for`), `:329-349`
  (`append_extension` only).
- **What:** only `<file>.bam.csi`, `<file>.bam.bai` and `<file>.cram.crai` are looked for. Picard
  and GATK write `sample.bai`. Without `--build-map-file-index` the run refuses with "missing
  index" although one is beside the file; with the flag it spends a full pass building a redundant
  `.csi` — which then sets up B5.
- **Proposal:** also try the replaced-extension name, as htslib does.
- **Downside:** none worth naming.
- **Effect:** none on calls; usability and speed. **Cost:** S. **Confidence:** high.

### B5 — A stale index is used without warning, and `.csi` is always preferred over a newer `.bai`

- **Where:** `src/bam/index_preflight.rs:306-317`.
- **What:** the first existing of `.csi`, `.bai` wins; modification times are never compared with
  the data file. **Scenario:** pvc once built `sample.bam.csi`; the user regenerates `sample.bam`
  and runs `samtools index` (writes `.bai`); pvc keeps reading the old `.csi`. Usually the offsets
  land mid-record and decoding fails loudly; if the old chunks end before the new data does, reads
  at a contig's tail are dropped silently (suspected). A stale `.crai` has the same exposure.
- **Proposal:** warn when an index is older than its data file; when both `.csi` and `.bai`
  exist, take the newer.
- **Downside:** copy tools that reset mtimes cause false warnings — hence warn, not refuse.
- **Effect:** output possible. **Cost:** S. **Confidence:** medium (selection traced; silent loss
  suspected).

### B6 — A leading soft clip lets an indel survive at the read edge that an unclipped read drops

- **Where:** `src/alignment/indel_norm.rs:335-350` (a clip gives the indel nowhere to shift, so it
  comes to rest right after the clip); the edge rule checks only CIGAR op index 0:
  `src/locus_generation/pileup/cigar_cursor.rs:339-340`, `src/locus_generation/pileup/decompose.rs:85-87`.
- **What:** `5S 3M 2D 40M` inside a homopolymer left-aligns to `5S 2D 43M`; the deletion is op 1,
  so it is kept, anchored on a base the read was clipped against. The same read unclipped
  (`3M 2D 40M` → `2D 43M`) has it dropped as op 0. The comments at `indel_norm.rs:154-157` and
  `:260-262` ("the cursor rejects first/last-op indels anyway") are false once a clip comes first.
  Effect: a deletion/insertion candidate at an un-normalised anchor inside the repeat, with support
  split from reads that place it correctly. Frequency unmeasured.
- **Proposal:** measure first (count prepared reads on HG002 chr20 whose first two ops are `S`
  then `I`/`D`). If material, decide "first/last" over non-clip ops in both the cursor and
  `decompose`, as the bad-CIGAR filter already does.
- **Downside:** changes output; the test `soft_clip_then_indel_away_from_chrom_start_keeps_indel`
  (`decompose.rs:426`) pins today's behaviour deliberately, so this is an owner decision.
- **Effect:** output yes. **Cost:** S (plus the measurement). **Confidence:** medium-high on the
  mechanism, low on the size.

### B7 — The kept-read window grows with the region length, and each read is deep-copied into it

- **Where:** `src/read/input/cursor.rs:685` (`kept.push_back(read.clone())`); eviction only in
  `move_to_region`, `:575-592`.
- **What:** every read emitted in a region stays held until the next move. The architecture doc
  (`doc/devel/ng/arch/alignment_cursor.md`, "Memory is bounded by depth plus the kept window")
  describes about 5 kb of overlap. For a generic segment of L bases at depth D the cursor holds
  about L·D/150 reads of ~500 B: ~100 MB for a 1 Mb segment at 30×, ~1 GB at 300×, per cursor.
  Each fresh read is also cloned once (four heap allocations) and again on every replay.
- **Proposal:** hold `Arc<AlignedRead>` in `kept`, and let the generator drop reads behind it
  mid-region (it already releases reference behind itself). First confirm how long generic
  segments get, especially with `--generic-path-only`.
- **Downside:** consumer signatures change (`locus_generation/ssr.rs`, the generator take owned
  reads).
- **Effect:** output none; memory and speed possibly better. **Cost:** M. **Confidence:** medium.

### Lower-consequence bugs

- **B8 — Prefetch thread panic leaves waiters blocked for ever.** `aligned_reads_reader/prefetch.rs:140`,
  `container.rs:512-523`: the claim stays `Decoding`. Fix: a guard that settles the claim on
  unwind. Turns a hang into an error. S, medium.
- **B9 — `MergedCursors::move_to_region` is not all-or-nothing as documented**
  (`src/read/input/sample_cursor.rs:193-196`). All callers treat the error as fatal, so only the
  doc is wrong today. S.
- **B10 — An interrupted index build leaves a truncated index at its final path**
  (`src/bam/index_preflight.rs:351-385`). Next run uses it and fails loudly; the user must delete it
  by hand. Write to a temporary name and rename. S.
- **B11 — `read.pos as u32`** (`src/read/prepared_read.rs:68`) wraps silently on a CRAM contig
  over 4.29 Gb. Use `u32::try_from`. S.

---

## Architecture

### A1 — Delete the dead reader stack in `src/bam`

- **Where:** `src/bam/segment_reader.rs` (1,896 lines) and `src/bam/segment_merge.rs` (551), both
  carrying `#![cfg_attr(not(test), allow(dead_code))]` (`segment_reader.rs:59`,
  `segment_merge.rs:28`) with "remove when the consumer lands"; `bam_input.rs`, `cram_input.rs`,
  `record_specs.rs` (a 6-line doc-only module); in `alignment_input.rs` `FilterCounts` (:150),
  `AlignmentMergedReaderConfig` (:251), `ContigInterval` (:342), `classify_pre_decode` and its
  enums (:368-433), `record_buf_to_mapped_read` (:442), `build_fasta_repository` (:846); the
  `AlignmentInputError` enum in `errors.rs`. Nothing outside these files references them.
- **What:** the consumer never landed; `src/read/input` was built instead. Each dead piece is a
  copy of a live one: `segment_merge` vs `sample_cursor` (k-way merge, duplicate detection),
  `classify_pre_decode` vs `read/filtering.rs` `verdict_on_raw_read` (same six-step flag/MAPQ
  cascade), `FilterCounts` vs `ReadFilterCounts`, `record_buf_to_mapped_read` vs
  `read/aligned_read.rs` `decode_record` (:117-180, line for line). About 3,500 lines including
  their tests. Their comments also name deleted things (`run_pileup`, `load_pileup_inputs`,
  `ia/feature_implementation_plans/…`).
- **Proposal:** delete them. Move `record_buf_to_mapped_read` into the one parity test that uses it
  as a reference (`read/aligned_read.rs`, around line 420), or drop that test. Repoint
  `doc/devel/ng/arch/alignment_file.md:295-310`.
- **Downside:** loses a reference implementation someone might want to compare against; git keeps it.
- **Effect:** none. **Cost:** M. **Confidence:** high.

### A2 — After A1, fold the rest of `src/bam` into `src/read`, and retire `MappedRead`

- **Where:** live items: `CigarOp`, `FLAG_*`/`DEFAULT_*` constants, `cigar_ref_span`,
  `cigar_is_bad`, `cigar_to_ops`, `read_exceeds_mismatch_fraction`, `compute_adaptor_boundary`
  (`alignment_input.rs`); all of `index_preflight.rs`; `AlignmentIndexError`; test-only
  `cram_files.rs`. Every production read goes `AlignedRead` → `into_mapped_read()`
  (`read/aligned_read.rs:281`) → `prepare_passthrough(MappedRead)` (`read/prepared_read.rs:62`),
  which sets a meaningless `source_file_index: 0` and copies `qual` (`:67`) where it could move it.
- **Proposal:** `prepare_passthrough(read: AlignedRead, …)` moving `qual`; delete `MappedRead`
  and `into_mapped_read` (4 call sites). Then move CIGAR/flag vocabulary to `read/cigar.rs` and
  `read/filtering.rs`, and index preflight to `read/input/index.rs`.
- **Downside:** ~20 `use` paths and doc links churn; `CigarOp` is also used by `alignment` and
  `locus_generation`, so `read/cigar.rs` becomes their dependency (it already is, via `bam`). The
  test-only `left_align_parity.rs` fixtures build `MappedRead` and need porting. Stopping after
  the `MappedRead` step is reasonable.
- **Effect:** none; one fewer copy per read. **Cost:** S for `MappedRead`, M for the move.
  **Confidence:** high.

### A3 — `read` → `psp`: move the read-filter provenance function to `run`

- **Where:** all four uses are in `src/read/filtering.rs`: `:124` (return type
  `Vec<(String, crate::psp::ParameterValue)>`), `:125` (import), `:174`
  (`crate::psp::header::wire_float_of`, a `pub(crate)` helper at `src/psp/header.rs:1142`), `:484`
  (test import). There is also a doc link to `crate::run::StoredSample` at `:115`.
- **What:** `ReadFilterConfig::provenance_parameters` turns filter settings into psp-header
  key/value pairs. Its only production caller is `src/run/gatherer.rs:661`; the keys are read back
  by prefix in `src/run/psp_caller.rs:48,1049`. `psp` has no production dependency on `read`.
- **Which way:** the delegated reviewer proposed moving it into `psp`. I recommend **`run`**
  instead (e.g. a `run/read_filter_provenance.rs`, or beside `gatherer.rs`): both the writer and
  the reader of these keys live in `run`, which already depends on `read` and `psp`. Moving it into
  `psp` only flips the edge — the file-format module would then know the read filter's settings.
  In `run`, `read` loses its `psp` dependency and `psp` gains none. `READ_FILTER_PROVENANCE_PREFIX`
  and `OFF` move with it.
- **Downside:** the exhaustive destructuring that makes a new filter field a compile error moves
  away from the struct; it still fails to compile because the fields are `pub`. `wire_float_of`
  must stay visible to `run`.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### A4 — `vcf` ↔ `read_phase` dependency cycle

- **Where:** `src/read_phase/stage.rs:36-37` imports `crate::vcf::{RecordPlace, line_fields}`;
  `vcf/writer.rs`, `vcf/mod.rs`, `vcf/assemble.rs`, `vcf/header.rs` import `read_phase`.
- **What:** `stage.rs` is a VCF-writer stage that parses and rewrites VCF text lines.
- **Proposal:** move `stage.rs` to `vcf/` (e.g. `vcf/read_phase_stage.rs`); `read_phase` then holds
  only data types and link scoring, with no outgoing dependency on `vcf`.
- **Downside:** the read-phase spec section no longer maps onto one module.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### A5 — psp: encoder and decoder share one description of the format (answer to the brief's question) — leave it

- **Where:** `src/psp/record.rs:185-206` (`RECORD_HEAD_FIELDS`), `:260-378` (`BODY_FIELDS` and
  the per-field name constants), `:465-479` (`record_fields`), `:521-612`
  (`RecordLayout::from_manifest`, `check_this_build_writes`), `:2867`
  (`the_fixture_encodes_to_these_exact_bytes`).
- **What I found:** every section's encoder and decoder live in the same file — `record.rs`
  (`encode_record_body_reusing` / `decode_record_body`), `footer.rs` (`encode_footer` /
  `decode_footer`), `index.rs` (`encode_index` / `decode_index`), `chain_ids.rs`
  (`LiveSetWriter` / `LiveSetReader`), `block.rs` (`BlockBuilder` / `BlockStream`). The field
  table is shared by both sides: the writer declares it into each file's header manifest, the
  reader refuses any file that declares a different list, and both sides name fields through the
  same constants. The table does *not* drive the codec — field order is written by hand twice —
  and the code says so (`record.rs:246-251`), naming the golden-bytes test that holds them in step.
  The writer's `finish` also decodes its own index and footer back before sealing the file
  (`writer.rs:782-822`), so the writer is checked by the reader's code.
- **Why leave it:** a table-driven codec cannot express which fields repeat per observation, and
  the record decoder runs at about twenty million records a second; an interpreter there would
  cost speed for no new safety. There is no writer/reader duplication worth merging in psp.
- **Small duplication I would also leave:** three byte cursors map the varint codec's two errors
  into their own error types (`block.rs:209-222`, `index.rs:216-243`, `record.rs:1595-1625`), about
  30 lines each. A generic cursor over an error trait would save ~60 lines and add a trait on the
  hot path. Not worth it.

### A6 — Production code kept alive only by tests in `psp`

- **Where:** `PspWriter::append` and its two helpers are `#[cfg(test)]` (`src/psp/writer.rs:409-583`),
  yet `append` is a documented spec feature with ~440 lines of tests (`writer.rs:2225-2662`) and
  three more uses in `src/run/psp_writer_line.rs:345,376,433` (tests). Likewise
  `SelectiveRecordIter` (`walk.rs:379-475`), `BlockRecords` (`block.rs:280-327`),
  `decode_record`/`encode_record_body` re-exports (`mod.rs:112-113`), and two header-only openers
  (`mod.rs:123,144`). 73 `#[cfg(test)]` attributes across `src/psp`.
- **What:** a reader cannot tell from the code which parts of the format API the caller uses; the
  append path in particular looks like a production feature with error variants
  (`PspWriteError::Reopen`, `UnsupportedHeader`) that no production run can raise.
- **Proposal:** decide per item. Thin test conveniences that compose production pieces
  (`decode_record`, `BlockRecords::split`) — keep. `append` — either delete it with its tests and
  error variants (nothing ships it), or say in its doc that it is a spec feature not yet wired.
  Recommend deletion; spec §6.4 can record it as not implemented.
- **Downside:** deleting `append` throws away tested code that a future resumable run might want.
- **Effect:** none. **Cost:** S–M. **Confidence:** medium (it is the owner's call whether append is
  coming).

### Other dependency directions (clean)

`psp` imports `locus_generation` (its record type), `types`, `region_typing`, `fasta`,
`reference_info` and a few others, and no `run`/`cli` (one doc link,
`src/psp/chain_ids.rs:120`). `bam` imports only itself. `alignment` imports `bam` (`CigarOp`),
`types`, `float`. No module in scope imports `cli`.

---

## Duplication

- **D1 — Five layers around one left-alignment call.** `LeftAlignPreparer<R, N>` →
  `AlignmentNormalizer` trait → `StructuredLeftAligner` → `left_align_indels` →
  `left_align_cigar(…, read_start = 0, remove_deletions_at_ends = false)`. Production uses only the
  default normaliser; the alternative (`FixpointLeftAligner`) is test-only, and the two fixed
  arguments make `LeftAlignResult.leading_deletion_bases_removed` always 0. Proposal: have
  `src/read/left_align.rs:183` call `indel_norm::left_align_indels` directly; drop the type
  parameter, the trait, the `Alignment` struct and `src/alignment/left_align_structured.rs`.
  Downside: a future aligner that moves the read start would need `Alignment` back; one example
  and one bench name `LeftAlignPreparer<R>`. Effect none. M. High.
- **D2 — The two STR delimiters repeat their anchor bookkeeping.**
  `src/alignment/ssr_best_path_flat_gap.rs:707-827` and `src/alignment/ssr_unit_robust.rs:898-1008`
  (~100 lines: junction offsets, flank-evidence counting, `resolve_anchors`, edge fallback);
  `best_of`/`UNREACHABLE` also twice. A shared traceback helper is possible; keep the recurrences
  separate. The flat-gap path is checked byte for byte by `delimit_parity`, so the change is safe
  to verify. Effect none. M. Medium.
- **D3 — Test helpers redefined per file** in `src/read/input` (`reference_bases`, `region`,
  `read_at` in `cursor.rs:1002-1075`, `sample_cursor.rs:351-363`, `open_bam.rs:1858`,
  `region_raw_aligned_reads.rs:332-340`, and four more) although `test_fixtures.rs` exists. Low
  value. S.

(The `src/bam` ↔ `src/read` duplication is A1.)

---

## Readability

### R1 — psp documentation describes deleted code and renumbered versions

- **Where:**
  - `src/psp/mod.rs:10-45` compares every term against "production's `src/psp/`" and "production's
    `PspReader`", which were deleted with the old caller on 2026-09-13.
  - `src/psp/record.rs:936-976`: `encode_record_body`'s doc says "The chain ids are dropped …
    Milestone E", but it calls `encode_record_body_reusing`, which encodes them
    (`encode_read_list`, around `:1063`). `:257` and `:50` repeat the claim.
  - `src/psp/record.rs:225` ("twelve of these twenty-five" — there are 26), `:254` ("the head's
    seven fields … arrive with C2" — there are 8), `:485` ("a twenty-eight-field declaration" —
    there are 34).
  - `src/psp/record.rs:176-178`: "It costs nothing today because no psp exists; from Milestone F
    it costs a version."
  - Field comments cite "format 1.2", "2.0", "2.1" (`record.rs:231-233, 1093, 1112, 2856-2921`,
    `walk.rs:1450-1451`, `block.rs:3196, 4214`) while `FORMAT_VERSION` is `(1, 1)`
    (`header.rs:497`), because the numbering was reset on 2026-10-18 (`header.rs:486-492`). A
    reader sees "added in format 2.1" in a 1.1 file.
- **Proposal:** delete the production comparison; fix the counts and the chain-id claim; replace
  version tags on fields with dates or drop them (the header doc already holds the history).
- **Effect:** none. **Cost:** S. **Confidence:** high.

### R2 — Stale documentation in `alignment` about which delimiter ships

`src/alignment/ssr_best_path_flat_gap.rs:5-31` says this aligner is "not the recommended
delimiter" and points to `ssr_best_path_unit_slip`, which no longer exists; in fact it ships as the
STR-bundle hull delimiter (`locus_generation/ssr_bundle.rs:81`). `ssr_unit_robust.rs:7,35,48`
cite three deleted modules; `src/alignment/mod.rs:16-28` lists two unimplemented algorithms as
implemented. Someone deciding whether the hull delimiter may change is told it is a frozen oracle.
S. High.

### R3 — `MarginalAligner` trait and its tests test only stand-ins

`src/alignment/mod.rs:486-576` (test-only trait, no real implementation) and the seven tests at
`:940-1199` check dummy implementations. Delete (~350 lines). S. High.

### R4 — History narrated in doc comments (taste; lower priority)

In `src/read/input` and `src/psp`, doc comments narrate milestones, review findings and dated
measurements (e.g. `src/read/input/cursor.rs:758-797`: ~36 comment lines for 2 lines of code;
`bam.rs`/`cram.rs` repeat the "positions, never bounds" essay that `aligned_reads_reader/mod.rs`
states). Keep the contract and the *why*; move history to `doc/devel/reports/`. This would roughly
halve the line count of several files. M. Medium — partly taste.

### R5 — Small items

- `src/alignment/mod.rs:260-262`: `RepeatGeometry.motif` must start at the tract's first base (the
  phase arithmetic at `ssr_unit_robust.rs:623` assumes it); say so.
- `src/bam/index_preflight.rs:215-226`: the mixed BAM+CRAM check cannot fire (files are passed one
  at a time from `open_bam.rs:255`); drop it or fix the doc that says it is enforced.
- `src/bam/alignment_input.rs:24-25`: `#[cfg(test)]` written twice.

---

## Looked at, left alone

- psp field table vs hand-written codec: deliberate, documented, pinned by a golden-bytes test (A5).
- psp header version check: higher minor versions are accepted and unknown trailing body fields
  skipped by encoding (`header.rs:711-726`, `record.rs:521-560`); per-observation extensions need a
  version raise, which the doc states.
- psp low-quality-read count cast `count as u32` (`record.rs:1438`): bounded by `num_obs` just above.
- psp truncation: `every_truncation_of_a_finished_psp_is_refused_without_panicking`
  (`reader.rs:1151`) covers every cut point.
- CRAM reference MD5: the vendored noodles fork checks each slice.
- BAM index chunks: noodles-csi merges overlapping chunks, so no record is read twice.
- Contig name/length/order/digest mismatch BAM vs reference: checked at open (`open_bam.rs:232-243`).
- Read groups: one declared group ignores record tags; several make a missing/undeclared `RG`
  fatal; a repeated `@RG ID` across the run is refused.
- Secondary/supplementary/unmapped drops are unconditional by design; QC-fail/duplicate are settings.
- Duplicate reads across one sample's files: caught on (position, flag, name) ties (`sample_cursor.rs:255-272`).
- MAPQ 255 becomes 0 and is dropped as low MAPQ (`aligned_read.rs:151`): counted, matches GATK.
- Prefetch cache contents and read order do not depend on thread timing.
- Merging the BAM and CRAM reader structs: they share only a `done` latch and a method shape; not worth it.
- `read_phase`: mates sharing a chain id de-duplicated, disagreeing mates removed, read cap
  deterministic; VCF stage timing and contig flush correct.
- `alignment` emission/stutter numerics: Q0 floored, NaN rates sanitised, −∞ only for impossible
  slips.
- `src/read/left_align.rs` does not duplicate left-alignment; it fetches the reference window and
  wraps the call (that wrapping is D1).
- `compute_adaptor_boundary` `u32` casts overflow only past 4.29 Gb, the crate-wide position limit.
