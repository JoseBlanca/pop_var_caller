# Code Review: ng_psp_mode_a1
**Date:** 2026-09-03
**Reviewer:** rust-code-review skill (orchestrator, 9 category sub-agents in isolated worktrees)
**Scope:** uncommitted diff for step A1 of `run_driver_psp_mode.md` — the psp header gains the run's identity-check fields
**Status:** Approve-with-changes

---

### 1. Scope

- What was reviewed: the working-tree diff on branch `ng-psp-mode` (base 9aa05795), exported as
  a patch and applied in each sub-agent's isolated worktree.
- In-scope files: [src/ng/psp/segmentation_section.rs](../../../../src/ng/psp/segmentation_section.rs) (new),
  [src/ng/psp/header.rs](../../../../src/ng/psp/header.rs),
  [src/ng/psp/mod.rs](../../../../src/ng/psp/mod.rs),
  [src/ng/psp/writer.rs](../../../../src/ng/psp/writer.rs),
  [src/ng/psp/record.rs](../../../../src/ng/psp/record.rs),
  [src/regions.rs](../../../../src/regions.rs),
  [src/ng/region_typing/mod.rs](../../../../src/ng/region_typing/mod.rs),
  benches/ng_psp_perf.rs, examples/dhat_ng_psp.rs, examples/ng_psp_parity.rs.
- Out of scope: all unchanged files; PROJECT_STATUS.md; the impl report's prose (verified
  separately at step 8a).
- Categories dispatched: reliability, errors, naming, defaults, idiomatic, refactor_safety,
  module_structure, smells (always-on set + multi-file scope), extras (parser of untrusted
  file input producing a stable format). unsafe_concurrency skipped (no unsafe/threads in the
  diff); tooling skipped (Cargo.toml untouched).

### 2. Verdict

**Approve-with-changes.** No behavioural defect in the shipped paths — the mutation pass shows
the *code* is right and several of the *tests* cannot see it going wrong. One Blocker (a fixture
whose own doc claim is false, leaving two of the three compatibility operands unguarded), one
process Major (the diff edits frozen `src/regions.rs` against the 2026-07-16 ruling), and a
cluster of test-coverage Majors, all with cheap fixes.

### 3. Execution status

- Run by the orchestrator in the dev container before dispatch:
  `cargo fmt --check` (clean), `cargo clippy --all-targets --all-features -- -D warnings`
  (clean), `cargo test --all-targets --all-features` (all 16 test binaries pass; lib suite
  "6044 passed; 0 failed; 14 ignored").
- Sub-agents re-ran scoped suites in their own worktrees (`ng::psp` 410 passed; `regions::`
  64 passed) and ran 12 mutations (5 survived, 0 changed-no-behaviour; every survivor proven
  behaviour-changing by a probe) plus 13 hostile hand-written TOML bodies through
  `Header::decode`.
- Not run: `cargo doc --no-deps`, `cargo audit` (no dependency changes), `cargo mutants`
  (not installed in the container; mutations were hand-applied and content-verified reverts).
- Findings labeled "Needs verification": 0.

### 4. Open questions and assumptions

1. **Where does `SegmentationInputs` live once psp and run depend on each other?** (affects M6)
   The arch's own §4 sketch puts `SegmentationInputs` inside the psp header, so the psp→run
   import is the design's, not an accident — but Milestone B adds the reverse arrow and the
   module_structure reviewer proposes lifting the type to `src/ng/segmentation_inputs.rs` with
   the existing `ng::run` re-export kept. Owner's call at Checkpoint A; nothing blocks on it.
2. **Is "same version, new required section" acceptable given no psp predates A1?** (affects Mi9)
   The plan's premise ("no ng psp exists outside tests yet") makes the pre-A1-file case empty;
   the cost is only that a *scratch* file from last week refuses with a damage-shaped message.
   Assumed resolved by the plan; flagged so the decision is recorded, not accidental.
3. **Are two catalog-less psps meant to compare equal?** The examples' "no catalog" sentinel
   (empty contig table, zero digest) makes any two such files agree on catalog identity when
   §6.2's check lands. Real `generate-psps` runs always have a catalog, so this only affects
   synthetic files; flagged for the checkpoint.

### 5. Top 3 priorities

1. **B1** — make the shared fixture's claim true (catalog `built-under` and `scan` differ from
   defaults): two small edits upgrade every existing round-trip test into the guard the doc
   already promises.
2. **M5** — un-edit frozen `src/regions.rs`: move the constructor into ng's own
   `GenomeRegions`, which only ever uses production's parser through its public API.
3. **M1** — close the forged-line hole for catalog contig names and tool-version, the exact
   bug class the header's field-name rule documents having fixed once already.

### 6. Findings

#### Blocker

- **B1: src/ng/psp/segmentation_section.rs:568 — the fixture's catalog `built_under` and `scan`
  are pure defaults, so a default-substituting decode passes the whole suite**
  **Categories:** reliability (mutation-proven M-A/M-B: both survived 410/410), defaults,
  refactor_safety (convergent). **Confidence:** High.
  The fixture doc says "every value differs from its type's default"; for two of the three
  operands of the future §6.2 refusal that is false, and mutations that replace the decoded
  sub-sections with defaults survive every test. Fix: perturb `built_under` (e.g.
  `min_flank_bp`, `min_score`) and all three `scan` fields in
  `segmentation_inputs_for_tests`; every round-trip test inherits the discrimination.

#### Major

- **M1: src/ng/psp/segmentation_section.rs (check_catalog) — a catalog contig name or
  tool-version carrying a newline forges a key line in the header text**
  **Categories:** reliability (probe-proven: the forged body round-trips equal), errors
  (convergent). **Confidence:** High.
  `check_contigs` and the manifest field-name rule refuse exactly this for their strings; the
  new section's writer-controlled strings have no such rule. Fix: refuse
  whitespace/control characters in catalog contig names and `tool-version` in
  `check_catalog`, with rule-table rows.

- **M2: src/ng/psp/segmentation_section.rs:415 — the empty-analysed-regions refusal has no test
  on either side** (reliability, mutation M-C survived 410/410). **Confidence:** High.
  Fix: the proposed `check_segmentation_refuses_empty_analysed_regions` test plus a reader-side
  body with zero `[[segmentation.analysed-region]]` rows.

- **M3: src/ng/psp/segmentation_section.rs:274 — the u64→u32 contig-bound clamp is untested;
  truncation survives because `i64::MAX as u32 == u32::MAX` coincides with the clamp**
  (reliability, mutation M-D survived 410/410). **Confidence:** High.
  At a contig length of 2^32 the truncating variant refuses every span the writer accepted.
  Fix: round-trip a span on a contig of length `1 << 32`.

- **M4: src/ng/psp/segmentation_section.rs:433 — the whole-contig analysed span (every
  whole-genome run's shape) never passes through the header in any test, and no segmentation
  value is pinned at exactly `MAX_TOML_INTEGER`**
  **Categories:** extras (probe-confirmed both boundaries currently accept), reliability
  (mutation M-E `>`→`>=` survived 410/410, convergent). **Confidence:** High.
  Every test fixture ends spans 50 bases short; the whole-contig fixtures live in
  benches/examples that never run under `cargo test`. Fix: a whole-contig round-trip test and
  a `min_flank_bp = MAX_TOML_INTEGER` acceptance (+1 refusal) test.

- **M5: src/regions.rs:171 — the diff edits frozen production code**
  **Categories:** module_structure. **Confidence:** High (the orchestrator re-read the ruling:
  [arch/typed_regions.md](../../ng/arch/typed_regions.md) Revision 2026-07-16 — "production
  (`src/ssr/`, `src/regions.rs`) is frozen; ng copies what it needs" — and
  [region_typing/mod.rs](../../../../src/ng/region_typing/mod.rs) itself still asserts
  "`src/regions.rs` is not edited (spec Revision)" a few lines above the diff's wrapper).
  Fix chosen (honours the freeze, no design change): revert `src/regions.rs` byte-identical to
  production; `GenomeRegions` stores its own span list, built through `RegionSet`'s public
  constructors on the BED/whole-genome paths and through ng's own validated constructor on the
  psp-decode path. Consumers use only `whole_contigs`/`from_bed_path`/`iter`/`len`/`is_empty`,
  so nothing else moves.

- **M6: src/ng/psp/header.rs:25 — the psp store now imports a run-stage type, and Milestone B
  adds the reverse arrow** (module_structure). **Confidence:** Medium — **deferred to the
  owner** (open question 1): the arch's §4 sketch itself places `SegmentationInputs` in the psp
  header, so this is the design's dependency; relocating the type's file is an arch edit, not
  implementer latitude.

#### Minor

- **Mi1: check_criteria accepts criteria `classify` release-asserts are impossible** —
  `min-purity = 7.0`, `bundle-threshold-bp = 0`, `period-max = 9` all decode believed (extras
  probes; idiomatic, errors, reliability convergent). Fix: enforce `classify`'s own three
  bounds in `check_criteria` (subsumes the NaN arm for purity).
- **Mi2: wire_float_of's `.expect` justification cites an invariant that does not hold on the
  `smuggle` path** (errors). The panic is impossible for a different reason — `f64::from_str`
  accepts everything `f32`'s `Display` produces, `NaN`/`inf` included. Fix the comment.
- **Mi3: the period-range refusal tests only `MinExceedsMax`, asserts only `is_err()`** —
  `ZeroMin` unreached (reliability). Fix: add `period-min = 0` and assert the field name.
- **Mi4: most new wire keys are unpinned — a serde rename would round-trip green and orphan
  every file** (extras). Fix: extend the body-spelling pin to every key the section writes.
- **Mi5: the ~18-line "no catalog" `SegmentationInputs` literal is pasted into the bench and
  both examples, and its zero-digest sentinel convention is documented nowhere**
  (smells, defaults, module_structure, idiomatic — 4-way convergent). Fix:
  `RepeatCatalogHeader::no_catalog(tool_version)` documented constructor; use it at all three.
- **Mi6: the `ContigIdentity → ContigBounds` clamp expression is spelled at four sites, its
  justifying comment on one** (smells, reliability, idiomatic, naming — convergent). Fix: one
  named helper for the in-crate sites.
- **Mi7: `for_period(u8::MAX)` probes `MinCopies`' private fallback via a sentinel that is only
  correct while `MAX_MOTIF_LEN < 255`** (smells, idiomatic). Fix: a `for_wider_periods()`
  accessor on `MinCopies` (an ng-owned type).
- **Mi8: encode-side field access is non-exhaustive — a field added in A2–A4 only errors at the
  decode literal, where a default can be quietly filled** (refactor_safety, probed compiling).
  Fix: exhaustive destructures in `From<&Header> for WireHeader` and the section's `from_*`.
- **Mi9: a required section landed in format 1.0, so a pre-A1 1.0 scratch file refuses as
  damage** (extras, errors, refactor_safety — convergent). **Resolved by the plan's premise**
  (no psp outside tests predates A1; the A1 impl report records the reasoning); no code change.
- **Mi10: `check_catalog` indexes parallel lists by loop counter after a length check** —
  `zip` makes the pairing structural (refactor_safety, probed). Fix: zip.
- **Mi11: naming** — `from_genomic_order_spans` understates its contract (rename
  `from_normalized_spans`); `SpanSetError::UnknownContig.contigs` is a count (rename
  `contig_count`); test binding `one_contig` holds two contigs; the criteria field is
  "the walk's" here and "the reader's" in segments.rs (naming). Fix during the M5 move.
- **Mi12: the catalog-md5 refusal does not name which contig row broke** (errors, smells).
  Fix: put the contig's name in the field path handed to `digest_of`.
- **Mi13: the fixture builds `start: 100, end: length - 50` and panics on contigs under 151
  bases** (reliability cross-note). Fix: derive the span from the contig's own length.

#### Nits

Missing `# Errors` heading on the new constructor (the sibling `from_bed_reader` has one);
`as_span` is fallible under the cheap-`as_` prefix (rename `to_span`); `previous`/`ordered`
bindings and the `class` abbreviation; a redundant `.to_string()` into `impl Into<String>`;
inline fully-qualified paths in the bench/examples where grouped imports would do; the
`eprintln!` in the 30k-scaffold test (kept deliberately: it carries the measured headroom
number and the harness captures it on passing runs); test local `three`.

### 7. Out of scope observations

- `SegmentationInputs::first_difference` ([segments.rs:166](../../../../src/ng/run/segments.rs#L166))
  compares fields non-exhaustively — a field added to the type would silently skip the very
  cohort check this header feeds (refactor_safety). Candidate fix in Milestone E, whose tests
  own that seam.
- Decode resolves each analysed span's contig name by linear scan — O(regions × contigs), which
  at the 35,000-scaffold assemblies the header cap budgets for is ~10⁹ comparisons worst case
  (idiomatic). Revisit if open-time ever shows up in a measurement; not on any hot path today.

### 8. Missing tests to add now

All carried inside findings B1, M2, M3, M4, Mi1, Mi3 (each with proposed code in the
per-category files under `tmp/review_2026-09-03_psp-header-a1/`), plus:
- `from_normalized_spans_accepts_every_normalized_set` — property test: whatever the BED
  parser's normalization hands out, the constructor accepts and returns verbatim, so the two
  definitions of "normalized" cannot drift (reliability).

### 9. What's good

- The two-sided rule set held: every hostile probe (13 hand-written bodies) was refused naming
  its field, and no allocation happens from an unchecked length
  ([header.rs:593](../../../../src/ng/psp/header.rs#L593) precedes all parsing).
- `check_contigs` running before span resolution is pinned by an existing message assertion —
  deleting the early call was caught (mutation M-J), a refactor-safety property most splits
  don't get for free.
- Decode rebuilds through the run's own checked constructors, so a value the run cannot build
  cannot arrive from a file either ([segmentation_section.rs](../../../../src/ng/psp/segmentation_section.rs)).
- The refusal tests sit one step over their boundaries (touching spans = gap of exactly zero;
  `chr2`'s span ends exactly at its length), which is what killed 7 of the 12 mutations.
- The intent check came back exactly clean: the §6.1 fields A1 owes, none of the
  amended-away fields (record count, boundary digest, `writer_version`) present (extras).

### 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- Scoped: `./scripts/dev.sh cargo test --lib 'ng::psp'` (was 410),
  `./scripts/dev.sh cargo test --lib 'region_typing'`
- Audit trail: per-category findings in `tmp/review_2026-09-03_psp-header-a1/`.
