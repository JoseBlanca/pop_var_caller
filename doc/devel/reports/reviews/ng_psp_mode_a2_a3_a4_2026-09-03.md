# Code Review: ng_psp_mode_a2_a3_a4
**Date:** 2026-09-03
**Reviewer:** rust-code-review skill (orchestrator, 9 category sub-agents in isolated worktrees; all nine were interrupted mid-run by a session rate limit and resumed in place — worktrees and applied patches survived)
**Scope:** uncommitted diff for bundled steps A2+A3+A4 over commit 114efe24 — the header's remaining §6.1 fields
**Status:** Approve-with-changes

---

### 1. Scope

- The diff: `Header.read_groups` (A2), `Header.observation_reach_ceiling` (A3),
  `WriterProvenance::record_read_filters` (A4), with fixtures across four psp files, two
  examples and the bench.
- In-scope files: src/ng/psp/header.rs (the substance), src/ng/psp/{mod,record,writer}.rs,
  src/ng/psp/segmentation_section.rs (one visibility change), examples/dhat_ng_psp.rs,
  examples/ng_psp_parity.rs, benches/ng_psp_perf.rs.
- Categories dispatched: the six always-on, plus defaults (serde defaults, the "off"
  convention), module_structure (a new dependency edge), extras (parser of untrusted input,
  stable format, intent check). unsafe_concurrency and tooling skipped (no triggers).

### 2. Verdict

**Approve-with-changes.** The wiring is right (the extras intent check found §6.1 complete
minus exactly the dropped record count; the reliability agent confirmed the rule table is
genuinely two-sided and the parity destructure forces classification of future fields), but
one test is decorative where it claims to guard, one function breaks the very guard pattern
this branch established, and the new dependency edge points the wrong way.

### 3. Execution status

- Orchestrator, pre-dispatch: fmt clean; clippy `-D warnings` clean; full suite 16 binaries,
  lib "6052 passed; 0 failed; 14 ignored".
- Sub-agents: reliability ran 9 mutations (6 survived, 0 changed-no-behaviour, every survivor
  probe-proven); refactor_safety probed its fix compiling and green; extras ran 7 hostile
  TOML bodies through `Header::decode` (all refused naming their field; the duplicated
  `@RG ID` accepted as intended).
- Findings labeled "Needs verification": 0.

### 4. Open questions and assumptions

1. **Which value does the real walk record as the reach ceiling?** The fixtures record the
   generator's *ceiling* (65,535); a real B1 walk will record its *configured*
   `max_record_span` (default 5,000). The fixtures' choice is honest for synthetic files
   whose records it truly bounds; B1's review should confirm the gatherer records the
   configured value. (naming cross-note)
2. Carried from A1 (unchanged): the format-1.0-unbumped question rides to Checkpoint A.

### 5. Top 3 priorities

1. **B1** — make the read-filter provenance test pin all six values with a config whose two
   booleans differ; today a qc-fail/duplicates transposition is invisible to the whole suite.
2. **M1+M2 together** — move the filter enumeration to `ng::read` as a method on
   `ReadFilterConfig` (fixing the psp→stage edge) and destructure the config exhaustively
   there (fixing the silent-unrecorded-filter hazard); psp keeps a generic parameter-recording
   seam.
3. **M3** — pin the deliberate acceptance of duplicated `@RG ID`s before someone "fixes" it
   into refusing real multi-file archives.

### 6. Findings

#### Blocker

- **B1: src/ng/psp/header.rs:1865 — `read_filters_land_in_provenance_with_off_spelled_out`
  cannot fail on four of the six recorded values** (reliability; four surviving mutants:
  recorded length off by one, fraction zeroed, and the qc-fail/duplicates sources transposed —
  the last invisible to the entire suite because every fixture has both booleans true).
  **Confidence:** High. Fix: pin all six values against a config whose booleans differ.

#### Major

- **M1: src/ng/psp/header.rs:211 — `record_read_filters` reads the config field-by-field, so
  a seventh filter compiles and goes unrecorded**, breaking the census-digest rationale its
  own doc states. **Categories:** errors, refactor_safety (probed: the exhaustive destructure
  compiles and passes), module_structure (convergent). **Confidence:** High.
- **M2: src/ng/psp/header.rs:211 — the new `psp → ng::read` import is a peer-imports-stage
  back-reference reaching the public surface** (module_structure). No header *field* holds
  the config — only this conversion — so it moves freely: put the enumeration on
  `ReadFilterConfig` in `ng::read` (stage-imports-infrastructure is the sanctioned
  direction, as with `ref_seq`), psp keeps a generic `record_parameters`. **Confidence:**
  High on the edge, Medium on the home (the recommended home is applied; the alternative,
  `ng::run`, separates the key-spelling contract from the format tests that pin it).
- **M3: src/ng/psp/header.rs:905 — the documented "duplicated `@RG ID` is legal" invariant
  has no test**; a later uniqueness tightening mirroring `check_contigs` would refuse real
  multi-file samples with the suite green. **Categories:** reliability, extras (convergent,
  both supplied the test). **Confidence:** High.

#### Minor

- **Mi1: `Bp(65_535)` spelled as a literal at three fixture sites, with a comment that
  mislabels it** — "the generator's shipped span cap" is `DEFAULT_MAX_RECORD_SPAN = 5_000`;
  65,535 is `MAX_RECORD_SPAN_CEILING`, the widest the generator *accepts*.
  **Categories:** defaults, naming, idiomatic, smells, refactor_safety — five-way convergent.
  Fix: the named constant, and the comment corrected.
- **Mi2: the walk-local number is `identifier` in the type, `walk-local-id` on the wire, and
  sits beside near-synonym `id`** (naming). Fix: rename the field `walk_local_id`.
  (idiomatic's alternative — drop the in-memory field and derive from position — was
  considered and not taken: the row is self-describing when it travels alone into E2's
  merge, and the both-sides redundancy check is this format's established style.)
- **Mi3: re-recording filters with the mismatch filter now off leaves a stale
  `read-filter-mismatch-bq-floor` beside `"off"`** (errors) — the misreading the conditional
  key exists to prevent; no such caller exists yet. Fix with M2: the read-side method
  documents the contract and names every key it owns so a re-recorder can clear first.
- **Mi4: the out-of-order fixture (`ReadGroupId(7)`) is far from the boundary** — a `!=`→`>`
  mutant survives, silently accepting entry 1 repeating identifier 0, the realistic
  two-pasted-walks shape (extras, probe-proven). Fix: add the identifier-0 case, keep the 7.
- **Mi5: the control-character rule is pinned only at `'\n'`** — narrowing to newline-only
  survived (reliability). Fix: a tab case.
- **Mi6: no reach-ceiling fixture at exactly `MAX_TOML_INTEGER`** — `>`→`>=` survived
  (reliability). Fix: extend the widest-number test.
- **Mi7: the `min_read_length > i64::MAX` digits-as-string fallback is untested** (errors,
  reliability convergent); a wrapping-cast regression would record a negative length. Fix:
  one assertion.
- **Mi8: `wire_float_of` widened to `pub(crate)` in the segmentation section while its
  family (`hex_of`, `digest_of`) lives in header.rs** (module_structure). Fix: move it home.
- **Mi9: two sentences in the step's impl report overstate what the tests reach**
  (reliability, step-8a class): "all three fields through every fixture" (the recorder has no
  caller outside header.rs) and "pins the on-values" (it pinned two of six). Fix the prose
  with the B1 fix.

#### Nits

Test locals `on`/`off`; the `record` closure reusing psp's most loaded noun as a verb;
`check_no_control_characters(field, value)` vs neighbour's `spelled`; `observation_reach_ceiling`
dropping the `_bp` suffix its sibling keeps (applied — cheap and mechanical); redundant
`.to_string()` into `impl Into<String>`; `"off"` spelled three times; inline qualified path in a
signature; whitespace-only `@RG ID` passing where `sample` would be refused (documented
asymmetry — verbatim SAM); one doc line noting the unconditional drops (secondary,
supplementary, unmapped) are not settings and so not among the recorded keys; plan checkboxes
A2–A4 flipped after commit per the loop's convention.

### 7. Out of scope observations

None new (A1's deferred items stand).

### 8. Missing tests to add now

Carried inside B1, M3, Mi4–Mi7; the reliability file supplies six as code.

### 9. What's good

- The rule table is genuinely two-sided: deleting the identifier-order check was caught in
  isolation with the exact expected message (reliability).
- The two-row read-group fixture catches an id/library column swap outright.
- The A1 exhaustive-destructure guards were extended everywhere for both new fields
  (refactor_safety verified all seven sites).
- `record_read_filters`' "off is a value, absence is a different fact" convention, with the
  one documented exception, is the defaults checklist done right (defaults).
- Hostile decode held: 7 hand-written bodies all refused naming their field (extras).

### 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib 'ng::psp'` (was 417) · full gate: fmt, clippy
  `-D warnings`, `cargo test --all-targets --all-features`.
- Audit trail: `tmp/review_2026-09-03_psp-header-a2a4/`.
