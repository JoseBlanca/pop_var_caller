# Fix Application Report: ng_read_filtering_stages_b1_2026-08-03.md

**Date:** 2026-08-03
**Source review:** `doc/devel/reports/reviews/ng_read_filtering_stages_b1_2026-08-03.md`
**Source state reviewed against:** branch `ng-generic-perf`, base `bfb54dd`, B1 working-tree diff
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
Blockers: 1 · Majors: 3 · Minors: 8 · Nits: 0

### Outcome totals
Applied: 12 · Applied with adaptation: 0 · Already fixed: 0 · Deferred: 0 · Disputed: 0
Failed validation: 0 · Blocked by context mismatch: 0 · Superseded: 0 · Awaiting user answer: 0

**Everything was applied.** Two items are recorded as deviations for the owner rather than
deferred, because the code change is in and it is the *design authority* that needs
reconciling — see §5.

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, **2,841 passed / 0 failed / 5 ignored**
- `cargo test --lib ng::` → 0, **1,542 passed / 0 failed / 2 ignored**
- `cargo test --examples` → 0, 52 passed / 0 failed
- `cargo doc --no-deps --lib` → 12 unresolved links, **all pre-existing**, none in a touched file
  (was 13 mid-run; the two links the review found were fixed)
- `cargo test --all-targets --all-features` → not run; pre-existing panic in
  `benches/psp_writer_perf.rs:386`
- Performance check → see §9; a walk-probe A/B was re-run because the fix adds one `open(2)`

**Four dumps byte-identical** to the `8cf6f03` baseline by `cmp`; walk probe anchor exact.

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | order-sensitivity untested | Apply | Applied | `open_bam.rs` | Pass + mutation |
| M1 | Major | equality ≠ resolvability; fail-fast lost | Apply | Applied | `open_bam.rs`, `mod.rs` | Pass + mutation + measurement |
| M2 | Major | one variant for two checks; message false | Apply | Applied | `open_bam.rs`, `mod.rs` | Pass |
| M3 | Major | operand order + prefix unasserted | Apply | Applied | `open_bam.rs` | Pass + mutation |
| Mi1 | Minor | `Reference` unproduced, doc false | Apply | Applied | `mod.rs` | Pass |
| Mi2 | Minor | three call sites bypass the check | Apply | Applied | `open_bam.rs` | Pass |
| Mi3 | Minor | count-mismatch branch unreached | Apply | Applied | `open_bam.rs` | Pass |
| Mi4 | Minor | `ReadFilterError::Reference` doc backwards | Apply | Applied | `filtering.rs` | Pass |
| Mi5 | Minor | two broken intra-doc links | Apply | Applied | `filtering.rs`, `cursor.rs` | `cargo doc` |
| Mi6 | Minor | self-contradicting module comment | Apply | Applied | `filtering.rs` | Pass |
| Mi7 | Minor | `SampleReads::cursor` doc silent on new failure | Apply | Applied | `mod.rs` | Pass |
| Mi8 | Minor | four impl-report defects | Apply | Applied | the impl report | N/A |

## 3. Questions asked and answers

None asked mid-run. Two items go to the owner at Checkpoint B as **deviations already applied**,
not as blocking questions — see §5.

## 4. Per-finding log

### B1 — the order-sensitivity had no test
- **Final status:** Applied
- **Reasoning:** the single most valuable finding of this review. Order-sensitivity is the
  *reason* the check is an equality rather than a resolvability test — both `over_records`' doc
  and `with_validated_contigs`' doc say so — and an order-insensitive rewrite passed all 1,538
  tests. On a perf branch, rewriting a 2,580-entry ordered walk as a hash lookup is a natural
  move, and it would have produced wrong variant calls with a green suite.
- **Implementation:** added `a_cursor_refuses_an_accessor_whose_contig_table_is_a_permutation`.
  It relies on `FIXTURE_CONTIGS`' two contigs having deliberately different lengths — its own
  doc says that is why.
- **Verification:** I re-ran the mutation in a **message-preserving** form (sort both tables by
  name, keep `first_disagreement`, so every rendered string is identical and only order-blindness
  changes), `grep`-confirmed present: `FAILED. 43 passed; 1 failed` — **the permutation test
  alone**. My first attempt used a hash-lookup mutation that also changed the message format and
  killed three tests; that would not have isolated order, so it was redone.
- **Residual risk:** None.

### M1 — table equality does not imply resolvability
- **Final status:** Applied
- **Reasoning:** the reviewer measured it, and it invalidated my "proves strictly more" claim.
  `ResidentRefSeq` and `WindowedRefSeq` take their `ContigList` as a constructor argument
  independent of the bytes, so a matching table can front a FASTA that cannot serve the contig.
  The deleted loop caught that; the first B1 did not.
- **Implementation:** a **third** check — one zero-length fetch on this cursor's own contig,
  after the table comparison. The loop asked that of every contig in the header for a property
  that matters only for the one about to be read; asking it once keeps the fail-fast at one
  `open(2)` instead of ~2,580. The check order is now argument → description → ability, which
  the probe also requires (it indexes `contig`).
  This revives `AlignmentFileError::Reference`, resolving Mi1 in the same stroke.
- **Verification:** added `a_cursor_refuses_an_accessor_that_cannot_serve_its_contig`, which
  drives a real `WindowedRefSeq` over a FASTA holding `chr1` only behind a table naming both,
  and asserts *both* directions — `chr1` accepted, `chr2` refused. Deleting the probe
  (`grep`-confirmed): `FAILED. 43 passed; 1 failed` — **that test alone**.
- **Residual risk:** None at the cursor. The constructors can still build a lying accessor;
  filed as an out-of-scope observation.

### M2 — one error variant for two checks
- **Final status:** Applied · **See §5, deviation 1.**
- **Reasoning:** the reviewer triggered the path and captured the message. Its headline —
  *"alignment file '…' does not match the reference contig table"* — is **false** once `open`
  has passed; what failed is the caller's accessor. Sharing one variant also left the two
  discriminable only by substring, which the first version of the new test demonstrated by doing
  exactly that.
- **Implementation:** `AlignmentFileError::CursorAccessorContigTable`, whose message names the
  accessor and the file. Both variants' docs now say how they differ and who is at fault.
- **Residual risk:** it deviates from spec §1 and §9 Q2 — §5.

### M3 — operand order and prefix unasserted
- **Final status:** Applied
- **Implementation:** the two assertions now pin the rendered shape —
  `contains("('chr1' vs 'not_chr1')")` and `contains("contig 'chr1': 100 vs 101")` — so the
  file-left decision is protected rather than merely recorded. (The prefix half of the finding
  dissolved with M2: the variant now carries what the prefix used to.)
- **Verification:** covered by the same runs as B1 and M1.

### Mi1–Mi8
All applied; each is one edit and none changes behaviour.
- **Mi1** resolved by M1 — `Reference` is produced again; its doc now scopes the fail-fast to
  this cursor's contig and explains *why* it is reachable (the constructor-argument gap).
- **Mi2** added `the_fixture_accessors_carry_the_same_contig_table_as_the_fixture_files`, a
  standing guard for the three `over_records` call sites that bypass `cursor`. Its doc records
  that this was false of every one of them before B1.
- **Mi3** added `a_cursor_refuses_an_accessor_whose_contig_table_is_shorter_than_the_files`.
- **Mi4, Mi6** rewrote two comments that explained themselves through deleted code; **Mi6**'s
  had contradicted another comment added by the same diff.
- **Mi5** the two dangling `[`ReadFilter::new`]` links de-linked. `cargo doc` back to 12.
- **Mi7** `SampleReads::cursor`'s doc now states that every accessor its factory hands out is
  checked, and what that costs a caller.
- **Mi8** the impl report corrected in four places, the largest being the deleted-test
  accounting: `read_filter_new_rejects_a_contig_missing_from_the_reference` asserted
  `UnknownContig` at a *count* mismatch, so **neither** its resolvability half nor its input
  class had moved. Both have successors now (M1's and Mi3's tests).

## 5. Deviations recorded for the owner (applied, not deferred)

1. **A new error variant — and the design authority is silent, not opposed.** §9 Q2's
   *illustrative* snippet reuses `ContigReconcile`; B1 adds `CursorAccessorContigTable`.
   **Correction, from B2's review:** this report first said it went against spec §1's "adds no
   new error". There is no such sentence — §1 says "change the meaning of any error", which
   adding a variant does not do, and the only "No new error type" statement is arch §4, scoped
   to `ReadFilterError`, a different enum. Applied because the reused variant produced a **false
   message**, demonstrated verbatim. Owner's call at Checkpoint B, on a smaller question than
   first stated.
2. **Spec §9 Q2's cost arithmetic is wrong**, confirmed independently by two agents: its 52 µs
   is a per-*accessor-construction* cost (34 µs of it a table clone) being multiplied as a
   per-*contig-fetch* cost. Not edited — a design document, and the claim is an argument for a
   change now made and measured.

## 6. Disputed findings
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** yes — M1's fix adds one `open(2)` per cursor, so the Checkpoint B measurement
  was re-run rather than inherited.
- **Result:** six runs each, same machine and session. Before (`bfb54dd`) mean **1.861 s**;
  after, as shipped, mean **1.834 s** — ≈27 ms per cursor. The pre-review implementation (no
  probe) measured 1.844; the gap between the two "after" variants is within the between-session
  spread, which is what one extra `open(2)` should look like. Honest figure: **~20–27 ms**.
- **Outcome:** pass. Consistent rather than noise — the slowest *after* run beats the fastest
  *before*.

## 10. Commands run
`cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`,
`cargo test --lib ng::`, `cargo test --lib ng::read::input::open_bam`, `cargo test --examples`,
`cargo doc --no-deps --lib`, `cargo build --release --examples`, the four acceptance dumps and
the walk probe against the `8cf6f03` baseline, and four mutations each `grep -c`-confirmed
present before its run and absent after the revert.

## 11. Command results
- fmt clean; clippy clean; `--lib` 2,841 / 0; `ng::` 1,542 / 0; examples 52 / 0
- `cargo doc` 12 unresolved links, all pre-existing
- four dumps **byte-identical**; walk probe anchor exact
- mutations: delete-comparison → killed; names-only → killed; **sort-both-tables → killed by the
  permutation test alone**; **delete-probe → killed by the resolvability test alone**
- final tree: `grep -c MUTATION src/ng/read/input/open_bam.rs` → 0

## 12. Notes
- **This review changed the design, not just the code.** B1 shipped as one check and now ships
  as three. The reviewers were right that a table comparison and a fetch prove different things,
  and that the first version quietly traded one for the other.
- **The Blocker is the fourth "test that cannot fail" this branch has surfaced in this plan**,
  and the first where the surviving mutation was a *plausible future optimisation* rather than
  an artificial edit.
- Three review agents ran in isolated worktrees; per-category findings are kept in the
  gitignored `tmp/review_2026-08-03_ng-read-filtering-stages-b1/`.
