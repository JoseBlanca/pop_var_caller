# Code Review: ng_paralog_filter_b3

**Date:** 2026-09-06
**Reviewer:** rust-code-review skill (orchestrator), three sub-agents in isolated worktrees
**Scope:** step B3 of the hidden-duplication filter plan — the writer's line entry point and the two-column patch
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `e3e0b6e5` on branch `ng-paralog-filter` — one new module
  (`patch.rs` and its tests), plus `write_line` and a public `RecordPlace` in the VCF writer.
- **Reviewed against:** `e3e0b6e5`. Every sub-agent detached its own worktree to that commit and
  confirmed both branch-only files before starting.
- **In-scope files:**
  - [patch.rs](../../../../src/ng/run/paralog_filter/patch.rs)
  - [tests.rs](../../../../src/ng/run/paralog_filter/patch/tests.rs)
  - the changed parts of [writer.rs](../../../../src/ng/vcf/writer.rs) — `write_line`,
    `RecordPlace`, `write_record` re-expressed
  - [tests.rs](../../../../src/ng/vcf/writer/tests.rs) — the five added tests
  - [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs) — the declaration and re-export
  - [ng_paralog_filter_b3_2026-09-06.md](../implementations/ng_paralog_filter_b3_2026-09-06.md) —
    its quantitative and mechanism claims
- **Deliberately out of scope:** `spill.rs` and `spill_file.rs` (reviewed at B1 and B2);
  `src/ng/paralog/` (Milestone A); `src/paralog/`, `src/var_calling/`, `src/sample_summary/`
  (production, frozen); the `WindowCoverage` stand-in in `mod.rs` (scheduled for deletion at the
  rebase onto `ng-window-coverage`).
- **Categories dispatched:** reliability + extras (a serialiser whose whole contract is
  byte-identity, on a per-record path at up to several thousand samples a line); errors + naming +
  idiomatic + defaults + the diff's own numbers; refactor_safety + module_structure + smells.
  `tooling` skipped — `Cargo.toml` untouched. `unsafe_concurrency` skipped — no `unsafe`, no
  threads, no shared state.

**One earlier fan-out against this same commit was killed before it wrote its findings files.**
Its three worktrees survived; their probe code was salvaged to
[salvaged_from_dead_agents/](../../../../tmp/review_2026-09-06_ng_paralog_filter_b3/salvaged_from_dead_agents/)
and handed to this round's agents to verify rather than re-find. **All ten salvaged probes
reproduced.** The most valuable of them is a nine-shape `every_shape()` fixture that encodes real
`VcfRecord`s through `record_line` — the round trip the committed tests do not have.

### 2. Verdict

**Request-changes.** One Blocker, seven Major, fifteen Minor.

**The code is correct on every shape ng emits** — the salvaged fixture proves it over nine encoded
record shapes, and a three-thousand-sample cohort line splices correctly. **What is missing is the
test that says so**, and five of the report's own claims about this step are wrong.

The single theme running through the Majors: **the step's structural claim is narrower than its
prose.** `write_record` really does go through one ordering check — verified structurally, no path
to the sink skips it. But the check reads a `place` that nothing ties to the line it accompanies,
and the *position* in that place still comes from a second copy of the padding rule.

### 3. Execution status

Run by the orchestrator in the dev container, on the branch at `8b731d41`:

| command | result |
|---|---|
| `cargo test --all-features --lib --bins --tests` | lib: `ok. 6427 passed; 0 failed; 15 ignored` in 44.43s; one integration test failed, pre-existing on `main` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 3 errors, all `needless_lifetimes` in `cohort_merge/build.rs:820`, `:893` and `serial.rs:67`, pre-existing |
| `cargo fmt --check` | dirty on 9 unique files, none of them this commit's, pre-existing |

The integration failure is
`ng_calling_loop_calls_genotypes::a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`
(`tests/ng_calling_loop_calls_genotypes.rs:1235`, left `[0, 0]` right `[0, 1]`).

`--all-targets`, `cargo doc --no-deps` and `cargo audit` were not run, for B1's and B2's reasons:
`examples/ng_candidate_selection_probe.rs` has not compiled against the current `ClosedLocus`
since `a33ada0f`, so `--all-targets` runs no test at all; the doc build is red on other modules'
links; and the container has no offline advisory database.

**All four failures are `main`'s at `a33ada0f` and none belongs to this step.** The gate the plan
uses instead — no failure `main` does not already have, and the `fmt` and `clippy` failure sets no
larger than `main`'s — holds.

**Findings labelled "Needs verification": 0.** Every finding below was produced by mutation or by a
probe run in a reviewer's own worktree, with output quoted in the per-category files under
[tmp/review_2026-09-06_ng_paralog_filter_b3/](../../../../tmp/review_2026-09-06_ng_paralog_filter_b3/).

**Mutation totals across the three agents: 15 run, 4 survived, 1 of those changed no behaviour** —
so three behaviour-changing survivors (M5 twice, M6 once), each proven to change behaviour by a
fixture that fails under the mutant and passes clean.

### 4. Open questions and assumptions

1. **Does `patch.rs` stay in `run/paralog_filter/`?** The implementation report's deviation 3
   argues it belongs there because "its rules are the filter's". Inventoried, the module encodes
   five rules — tab separates columns; `FILTER` is seventh and `INFO` eighth; `.` is missing;
   `PASS` means every filter passed; `;` joins both columns — and **all five are VCF grammar**.
   The filter id and the two `INFO` keys arrive as arguments. Affects **Mi1**, **Mi2**. Reversing a
   recorded deviation from a committed step is not the implementation loop's call; raised at
   Checkpoint B, with the minimum fix (import the vocabulary rather than re-spell it) applied
   meanwhile.
2. **Should `write_line` refuse a line containing a newline, and an added value containing a tab?**
   Neither is reachable from `record_line`, and pass three's added values are C4 constants — but
   the spill reader validates a length cap and not the bytes, so a corrupt frame yields a "line"
   holding the next entry's binary payload. Affects **Mi10**. The alternative to refusing is a doc
   sentence saying which malformednesses are *not* refused.

### 5. Top 3 priorities

1. **B1** — the step's own contract has no committed test; all 13 fixtures are hand-written lines,
   and they are not the shape the encoder emits.
2. **M1 + M3** — `write_line`'s ordering guarantee rests on a `place` nothing ties to the line, and
   the field C4 must supply is a bare `u64` where the spill entry carries a `Position`. One fix
   closes both, and it is cheapest now: `RecordPlace` has 7 sites, all inside `writer.rs` and its
   tests.
3. **M2** — the padding rule has two copies, so `write_record`'s ordering key and its own `POS`
   column can drift. Measured: mutating one copy fails 1 of 20 writer tests.

### 6. Findings

#### Blocker

- **B1: [patch/tests.rs:1-205](../../../../src/ng/run/paralog_filter/patch/tests.rs#L1-L205) — the step's contract is tested only against lines the encoder does not emit**
- **Categories:** reliability (filed), errors and refactor_safety (cross-category notes) — convergent
- **Confidence:** High
- **Problem:** Every fixture is a hand-typed byte literal. `a_record_line()` (line 11) is
  `b"SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5;AC=1\tGT:GQ:AD\t0/1:30:5,5\t0/0:99:8,0"`, and
  three of its columns are shapes the encoder never writes: `AF=0.5` is two decimals where
  [encode.rs:36](../../../../src/ng/vcf/encode.rs#L36) fixes six; `GT:GQ:AD` is a `FORMAT` string
  [encode.rs:343](../../../../src/ng/vcf/encode.rs#L343) never produces (it writes `GT:GQ:DP:AD`);
  and `AF=0.5;AC=1` puts `AC` where a real line carries `AN` and `DP`. Nothing in the committed
  suite builds a `VcfRecord`, encodes it through `record_line`, and asserts the round trip.
  Spec §10's standing oracle — filter on at an unreachable target, the two `INFO` keys stripped,
  equals the filter-off file — is therefore checked against lines ng does not write.
- **Why it matters:** If the splice ever mis-handles a shape the encoder produces, the output is a
  VCF with shifted sample columns — wrong genotypes for every sample, silently, with no panic and
  no failing test. The rubric names absence of tests for exactly this class.
  **The behaviour is right today**: the salvaged `every_shape()` fixture runs nine encoded shapes
  through the patch and all nine come back byte for byte. What is missing is the test that would
  notice when a later change to `encode.rs` breaks it.
- **Suggested fix:** commit `every_shape()` and its round-trip assertions into `patch/tests.rs`,
  extended two ways. First with the three shapes the salvaged fixture misses (M-list item below):
  a record on a contig other than the first, a repeat tract whose sample is a no-call (`REPCN` is
  `.`), and a three-thousand-sample cohort. Second — and this is the half both salvaged probes
  skip — **assert what `FILTER` and `INFO` become**, not only that the other columns did not move:
  both probes contain `if index == 6 || index == 7 { continue; }`, so across all nine shapes no
  assertion constrains the two columns the function exists to change. An implementation that wrote
  the ratio into `FILTER` and the tag into `INFO` passes both. See *Missing tests* items 1 and 2,
  which were written and run green (`test result: ok. 29 passed; 0 failed`).

#### Major

- **M1: [writer.rs:144](../../../../src/ng/vcf/writer.rs#L144) — the ordering check reads a `place` nothing ties to the line, so `write_line` can write a backwards VCF and return `Ok`**
- **Categories:** reliability, errors, refactor_safety — convergent across all three
- **Confidence:** High
- **Problem:** `write_line(place, line)` calls `check_order(last, place)` and never reads the
  line's bytes. Reproduced independently by two agents: two calls with ascending places and
  descending `POS` columns both return `Ok`, and the finished file's positions come out
  `["900", "100"]`. The same holds for tract-ness — a place that lies about it admits a tie the
  writer exists to refuse. The method's own `# Errors` section
  ([writer.rs:140-143](../../../../src/ng/vcf/writer.rs#L140-L143)) says "If the **line** runs
  backwards", naming a subject the check never inspects, and the module doc
  ([writer.rs:4-6](../../../../src/ng/vcf/writer.rs#L4-L6)) says ordering "is checked here".
- **Why it matters:** spec §6 trap 6 is *bypassing `write_record` must not bypass its check*. The
  check is not bypassed — but it can be fed a place unrelated to the bytes, which produces the
  backwards VCF the trap names. The commit message's "There is nowhere left to bypass" is true of
  the check and false of the guarantee the check exists to give. Not reachable today (nothing calls
  `write_line`), which is what makes now the cheap moment.
- **Suggested fix:** make the mismatch hard to build rather than documenting it. With M3's change
  landed first so the types line up, add to `paralog_filter/mod.rs` (the stage importing the format
  module — the correct direction):

  ```rust
  impl From<&SpillEntry> for RecordPlace {
      fn from(entry: &SpillEntry) -> Self {
          Self {
              at: GenomePosition { contig: entry.contig, position: entry.position },
              is_repeat_tract: entry.is_repeat_tract,
          }
      }
  }
  ```

  C4 then writes `writer.write_line((&entry).into(), &patched)?` and never spells a field. Fix the
  `# Errors` wording to say *place* rather than *line*, and state what the writer cannot check.
  Add *Missing tests* item 5 as a pinning test either way.

- **M2: [writer.rs:231](../../../../src/ng/vcf/writer.rs#L231) and [encode.rs:154](../../../../src/ng/vcf/encode.rs#L154) — the padding rule has two copies, so a record's ordering key and its own `POS` column can drift**
- **Categories:** refactor_safety
- **Confidence:** High
- **Problem:** `place_of` and `written_position` are the same six lines, spelled twice:

  ```rust
  match record.padding_base() {
      Some(PaddingBase::Left(_)) => start - 1,
      Some(PaddingBase::Right(_)) | None => start,
  }
  ```

  `write_record` ([writer.rs:123-126](../../../../src/ng/vcf/writer.rs#L123-L126)) now calls
  **both** in one expression — `record_line` writes `POS` from `written_position`, `place_of`
  builds the ordering key — so the two halves of one record come from two independent spellings of
  one rule. Measured, by mutating `place_of`'s left-padding arm to `=> start` and leaving
  `encode`'s copy intact: `FAILED. 19 passed; 1 failed`, the single failure being
  `a_tract_followed_by_a_generic_locus_at_one_position_is_refused`. One test of twenty, against the
  nine of twenty the report quotes for deleting the ordering check. Notably
  `the_order_is_checked_on_the_written_position_not_the_span_start` did **not** fail — it pins
  `place_of` against the span start, not against the `POS` column the same record gets.
  **The duplication pre-dates this commit**; `place_of`'s body is untouched by `e3e0b6e5`, only its
  call site moved. What is this step's is the claim, in deviation 1 and in the commit message, that
  the re-expression gives "one ordering check rather than two that could drift" — true of the
  check, not of the position it reads.
- **Why it matters:** it is the hazard of M1 one level down, arriving not from a careless call site
  but from an edit to one of two copies, with one test in the way.
- **Suggested fix:** give the rule one owner — make `written_position` `pub(crate)` and have
  `place_of` call it. Verified by the agent: `cargo test --all-features --lib "ng::vcf"` →
  `ok. 137 passed; 0 failed`. `PaddingBase` then becomes unused in `writer.rs`'s import.

- **M3: [writer.rs:53](../../../../src/ng/vcf/writer.rs#L53) — `RecordPlace.position` is a bare `u64` where the spill entry beside it carries a `Position`, and the report claims they match**
- **Categories:** naming, refactor_safety — convergent
- **Confidence:** High
- **Problem:** Deviation 2 and the commit message say `RecordPlace` "is **exactly** the three head
  fields the spill entry carries". Two of the three match:

  | `RecordPlace` ([writer.rs:54-59](../../../../src/ng/vcf/writer.rs#L54-L59)) | `SpillEntry` ([spill.rs:139-145](../../../../src/ng/run/paralog_filter/spill.rs#L139-L145)) |
  |---|---|
  | `pub contig: ContigId` | `pub contig: ContigId` |
  | `pub position: u64` | `pub position: Position` |
  | `pub is_repeat_tract: bool` | `pub is_repeat_tract: bool` |

  So C4, holding a `SpillEntry`, must write `position: entry.position.get()` — an untyped unwrap at
  precisely the call site where M1's mismatch is accepted, and a field-by-field literal of two
  similar-looking numeric fields with nothing to check them against. A constructor does not fix it:
  `RecordPlace::new(contig, position, is_repeat_tract)` takes the same three values positionally.
  Separately, `check_order` hand-rolls `(next.contig.0, next.position) < (last.contig.0, last.position)`,
  reaching through `ContigId`'s tuple field, where
  [`GenomePosition`](../../../../src/ng/types.rs#L59-L63) derives `Ord` and its own doc says it
  exists to "serve directly as a sort key wherever reads or loci are ordered along the reference".
- **Why it matters:** this is the field pass three's ordering guarantee rests on, and the two types
  that are supposed to agree do not. `grep -rn RecordPlace src/ tests/ examples/` returns **7 sites,
  all in `writer.rs` and `writer/tests.rs`** — after C4 lands it returns more.
- **Suggested fix:** carry the crate's genome-order key. Verified by two agents independently
  (20 of 20 writer tests pass; with M1, M2 and Mi5 also applied, `cargo test --all-features --lib`
  → `ok. 6428 passed; 0 failed; 15 ignored`):

  ```rust
  pub struct RecordPlace {
      /// Which base the record is written on, after the padding rule has moved it.
      pub at: GenomePosition,
      /// Whether it is a repeat tract, which is what admits the one legal tie.
      pub is_repeat_tract: bool,
  }
  ```

  `check_order` then collapses to `if next.at < last.at` and `let tied = next.at == last.at`. One
  edit neither probe mentions: `writer/tests.rs` must import `ContigId` explicitly, because it
  currently inherits it through `use super::*` from a `use` line this change removes.

- **M4: [patch.rs:119-136](../../../../src/ng/run/paralog_filter/patch.rs#L119-L136) — the error's message, its doc, the `# Errors` section and the field doc each describe a different condition**
- **Categories:** errors, reliability — convergent
- **Confidence:** High
- **Problem:** Four statements about one condition, no two alike.
  - The **code** refuses when `splitn(9, …)` yields fewer than nine pieces, so a **nine-column**
    line — the eight fixed columns plus `FORMAT`, no samples — is accepted and patched. Reproduced:
    `chr1\t100\t.\tA\tG\t50.0\thiddenParalog\tAN=2;DP=10\tGT:GQ:DP:AD` comes back patched.
  - The **variant doc** says the opposite: "A record's line always has at least ten … Fewer means
    the line did not come from the encoder."
  - The **`# Errors` section** says a third thing: "If the line does not have the eight columns a
    VCF record has before `FORMAT`." An eight-column line *has* those eight and is still refused.
  - The **field doc** adds a fourth: `/// How many pieces the line split into, at most {PIECES}.`
    `{PIECES}` is not interpolated in a doc comment — rustdoc renders it literally — and the bound
    is wrong regardless: the variant is only built inside `pieces.len() < PIECES`, so `columns` is
    at most **8**.

  The invariant the doc states is real — `VcfRecord::new` asserts `!sample_columns.is_empty()`
  ([vcf/mod.rs:232-233](../../../../src/ng/vcf/mod.rs#L232-L233)) — which makes the doc right and
  the guard one short.
- **Why it matters:** C4 is the first caller, and a reader of the doc will believe a truncated
  nine-column line is refused. It is not, and `write_line` validates nothing either, so the two go
  out together. A nine-column record line in a cohort VCF declares samples in the header and gives
  none.
- **Suggested fix:** pick the rule and state it once. Either tighten the guard to `< PIECES + 1`
  and add `a_line_with_a_format_column_and_no_samples_is_refused`, or say nine everywhere and
  explain why the tenth is not checked. Fix the field doc to "at most eight" and drop the
  uninterpolated `{PIECES}`.

- **M5: [patch.rs:87](../../../../src/ng/run/paralog_filter/patch.rs#L87), [:96-99](../../../../src/ng/run/paralog_filter/patch.rs#L96-L99) — two branches have no test at all; both mutations survive all 13**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Dropping `&& !filter.is_empty()` from `write_filter` leaves all 13 patch tests green,
  and turns an empty `FILTER` column into `;hiddenParalog`. Deleting the `info_to_add.is_empty()`
  early return from `write_info` also leaves all 13 green, and makes a `.` `INFO` with nothing to
  add come back as an **empty column** — breaking the byte-identity the module doc
  ([patch.rs:47-49](../../../../src/ng/run/paralog_filter/patch.rs#L47-L49)) promises for that
  exact case.
- **Why it matters:** "nothing to add returns the line unchanged" is the property deviation 4 calls
  this step's own proof, and it is tested at exactly one `INFO` spelling.
- **Suggested fix:** *Missing tests* items 6 and 7, both written and run.

- **M6: [writer.rs:144-159](../../../../src/ng/vcf/writer.rs#L144-L159) — nothing says what the writer's state is after it refuses a line**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Advancing `self.last` and `self.records_written` *before* `check_order` survives all
  20 writer tests. Under that mutant a refused line is counted in `records_written()` and becomes
  the place the next line is checked against, so a line behind the last one actually written is
  then accepted. Every existing refusal test stops at the `Err`.
- **Why it matters:** `records_written()` is the run report's number, and the ordering check is only
  as good as the state it carries. Both would be wrong with no test failing.
- **Suggested fix:** *Missing tests* item 3, run against the mutant and failing with
  `a line the writer refused was counted as written; left: 2, right: 1`.

- **M7: [patch.rs:132-135](../../../../src/ng/run/paralog_filter/patch.rs#L132-L135) — the error names no record, so a pass-three failure cannot be traced**
- **Categories:** errors
- **Confidence:** High
- **Problem:** `TooFewColumns { columns }` carries a count and nothing else. Pass three walks the
  spill calling this once per record, so an operator whose run dies on a file of millions of records
  sees `a record's line splits into 3 column(s) …` with no contig, no position, no ordinal and no
  prefix of the offending bytes — all of which the spill entry beside it carries.
- **Why it matters:** this is the one error the filter's third pass can raise.
- **Suggested fix:** `rewrite_filter_and_info` genuinely does not know where the line came from, so
  the right home is C4's wrapper. Say so in the doc now, so C4's author knows it is their job:
  "**The error names no record** — this function is handed bytes and nothing else. A caller walking
  a spill should wrap it with the entry's contig and position, or the failure is untraceable."

#### Minor

- **Mi1: [patch.rs:1](../../../../src/ng/run/paralog_filter/patch.rs#L1) — all five rules the module encodes are VCF grammar, so deviation 3's argument for its location does not hold.** Inventoried: tab separates columns; `FILTER` is seventh and `INFO` eighth; `.` is missing; `PASS` means every filter passed; `;` joins both. None is the filter's — the id and the two keys are arguments. A future annotator (`--exclude-regions`, a post-hoc `FILTER`) would import from `run::paralog_filter` to reuse it. **Deferred to Checkpoint B** (open question 1); the minimum, Mi2, is applied meanwhile.
- **Mi2: [patch.rs:40,43](../../../../src/ng/run/paralog_filter/patch.rs#L40-L43) — `MISSING` and `PASS` re-spell values the crate already owns.** `crate::ng::vcf::MISSING_FIELD` ([encode.rs:21](../../../../src/ng/vcf/encode.rs#L21)) and `FilterVerdict::Pass.as_str()` ([vcf/mod.rs:706](../../../../src/ng/vcf/mod.rs#L706)) are the same two values, and the latter is what the header declares. Three copies of `"."` and two of `"PASS"` kept in step by nothing. Verified: taking both from the parent compiles once `as_str` is `const fn`, 13 of 13 pass.
- **Mi3: [vcf/mod.rs:54](../../../../src/ng/vcf/mod.rs#L54) — `RecordPlace` is in a `pub` method's signature but not re-exported**, while every other type in a `pub` signature in that module is. C4 would import it one level deeper than `VcfWriter` for no visible reason. Fix: `pub use writer::{RecordPlace, VcfWriteError, VcfWriter};`
- **Mi4: [writer.rs:163](../../../../src/ng/vcf/writer.rs#L163) — `check_order` reads `RecordPlace`'s fields one at a time**, so a fourth component of the ordering key would not surface in the one function whose correctness depends on the whole field set. `place_of` already builds the literal exhaustively, so construction is covered and consumption is not.
- **Mi5: [patch.rs:59-66](../../../../src/ng/run/paralog_filter/patch.rs#L59-L66) — length-check-then-index, where the guard and the indexing are coupled by hand.** Measured: an off-by-one in the guard is an index **panic**, not an error — `index out of bounds: the len is 8 but the index is 8` at `patch.rs:77`. **Take the array form, not the slice form.** Both salvaged probes use a slice pattern, whose arity is checked at run time: measured at `PIECES = 10` the slice form compiles and fails 10 of 13 tests at run time, while a fixed-size array via `try_into` gives `error[E0527]: pattern requires 9 elements but array has 10`. Only the array form makes the compiler flag the refactor.
- **Mi6: [patch.rs:54](../../../../src/ng/run/paralog_filter/patch.rs#L54) — the owned return is what makes step D2's deferred buffer question cost two changes instead of one.** Deferring the reuse decision is right; the signature can absorb only one answer. Adding a `rewrite_filter_and_info_into(&mut Vec<u8>, …)` form now, with the owning one as a wrapper, leaves D2's measurement free and changes no call site. Sized: the tests' own sample column is 11 bytes with its tab, so a 3,000-sample line is about 33 kB allocated and copied per record.
- **Mi7: [patch.rs:82,95](../../../../src/ng/run/paralog_filter/patch.rs#L82) — `write_filter` and `write_info` do not write.** They append to a `Vec<u8>`; their nearest neighbours in the crate — `VcfWriter::write_line`, `write_record`, `write_stream`, `Sink::write_all`, all in a file this same commit changes — write to a file, and one can fail doing it. `append_filter_column` / `append_info_column`, or `push_*`.
- **Mi8: [patch.rs:120-136](../../../../src/ng/run/paralog_filter/patch.rs#L120-L136) — `#[non_exhaustive]` is on the enum but not on its one struct variant**, so adding a field to `TooFewColumns` — which M7 asks for — is still a breaking change for a downstream `match`.
- **Mi9: [patch.rs:54](../../../../src/ng/run/paralog_filter/patch.rs#L54), [mod.rs:36](../../../../src/ng/run/paralog_filter/mod.rs#L36), [writer.rs:53,144](../../../../src/ng/vcf/writer.rs#L53) — four new items are `pub` to the world where `pub(crate)` is what the design needs.** All four have one intended caller inside the crate, and at this commit none at all. Counterpoint recorded: `spill` and `spill_file` beside it are `pub mod` too, and so is every module on the path, so narrowing `patch` alone makes it the odd one out — if the surface is worth narrowing, narrow the `paralog_filter` subtree as one decision at the end of the milestone.
- **Mi10: [patch.rs:54-79](../../../../src/ng/run/paralog_filter/patch.rs#L54-L79) and [writer.rs:144](../../../../src/ng/vcf/writer.rs#L144) — neither refuses a tab or a newline, and one of either silently changes the file's shape.** A `0x0A` in a line makes one `write_line` call put two record lines in the file while `records_written()` counts one, the second possibly behind the first. A `0x09` in an added value adds a column, shifting every sample column right. Neither is reachable from `record_line`, but the spill reader validates a length cap and not the bytes. The report's "The patch is what refuses a malformed one" (line 143) is broader than the code: it refuses one malformedness of three. See open question 2.
- **Mi11: [patch.rs:54](../../../../src/ng/run/paralog_filter/patch.rs#L54) — not idempotent, and no doc or test says once-only.** A second call gives `FILTER = hiddenParalog;hiddenParalog` and a repeated `PARALOG_LR`, which strict readers reject. `the_filter_and_the_info_are_independent` establishes that the two columns can be patched separately, which is the reading under which a caller makes two calls. Correct as is — the fix is one doc sentence and the pinning test.
- **Mi12: [patch/tests.rs:91-98,111-121](../../../../src/ng/run/paralog_filter/patch/tests.rs#L91-L98) — three of the 13 tests pin inputs the encoder cannot produce, and are presented as the tests of the two easy-to-reverse rules.** `FilterVerdict`'s five values are all words, so `FILTER` is never `.` and never empty; `info_column` pushes `AN=` and `DP=` unconditionally, so `INFO` is never `.` and never empty — confirmed over all nine encoded shapes, the thinnest being `AN=0;DP=0`. Keep the defensive behaviour; retitle to say it is defensive, and cover the *reachable* half of trap 7 (a real `notPeriodic` or `EMNoConv` line being joined to) on an encoded line.
- **Mi13: [patch.rs:112-117](../../../../src/ng/run/paralog_filter/patch.rs#L112-L117) — `added_bytes` has no test and its only effect is invisible to every test.** Replacing the body with `0` leaves all 29 patch tests green — the one mutation run that changed no behaviour, since it only sizes `Vec::with_capacity`. Probed over 36 shapes it is never an under-estimate and over-estimates by at most 5 bytes, but nothing states that it must be an upper bound on what `write_filter` and `write_info` write.
- **Mi14: [patch.rs:54](../../../../src/ng/run/paralog_filter/patch.rs#L54) — no property test on a function that is a parser, a serialiser, and an identity law.** The reliability checklist requires one for exactly this shape. `proptest` is already a dev-dependency and `proptest-regressions/` exists.
- **Mi15: [writer/tests.rs:408-427](../../../../src/ng/vcf/writer/tests.rs#L408-L427) — two of the five new writer tests cannot fail on the most obvious way to break `write_line`.** `a_line_written_without_its_record_reaches_the_file_exactly` asserts `ends_with`, which a writer that duplicates the line, or writes anything extra before it, or drops the header, still satisfies. `writing_a_record_and_writing_its_line_produce_the_same_bytes` compares two files both produced through `write_line`, so any defect in `write_line` cancels. Measured: a mutant writing every line and newline twice fails **none** of the five new tests — it is caught only by two written before this step. Fix: compare the whole body, not a suffix.

#### Nits

- `added_bytes` returns a count but reads as though it returns bytes; `added_byte_count` says which.
- `PIECES` names a shape, not a domain thing; `COLUMNS_THE_VERDICT_TOUCHES` is used once, only to
  compute `PIECES`, while `INFO` is reached with the unnamed `pieces[COLUMNS_BEFORE_FILTER + 1]`.
  With Mi5's array destructure both constants stop indexing anything — collapse to a single
  `PIECES` whose doc carries the sentence they were carrying.
- `patch.rs:5-6` says the sample columns are "all but **eight** of a cohort's line"; they are all
  but **nine** — the eight fixed columns plus `FORMAT`. `patch/tests.rs:5` gets it right.
- `write_info` branches on `says_something || index > 0` inside the loop, a loop-invariant plus an
  index test where a running `needs_separator` flag reads more directly.
- `RecordPlace` derives `Clone, Copy, PartialEq, Eq, Debug` but not `PartialOrd, Ord`, which is why
  `check_order` compares by hand; M3 gets them from `GenomePosition`.
- `VcfWriteError` is public and **not** `#[non_exhaustive]`, while `LinePatchError` (this commit)
  and `SpillError` both are. Pre-existing; this commit makes the inconsistency visible by adding
  the third.

**The 18 test names: no findings.** All read as sentences saying what must hold, none is a
placeholder, none abbreviates, and none needs the body read to know what it asserts.

### 6a. The diff's own numbers

Every claim re-derived by running something. **Eleven CHECKED-CORRECT, five WRONG.**

| claim | verdict | the real value |
|---|---|---|
| `patch.rs` 139 lines; `patch/tests.rs` 205 lines, 13 tests; `writer.rs` +35/−7 = 28 net; `writer/tests.rs` +116 | CHECKED-CORRECT | as stated (`git show --numstat`) |
| 13 patch tests; 5 new writer tests, 20 total (15 before) | CHECKED-CORRECT | as stated |
| mutation 1 → **9 of 20** writer tests fail | CHECKED-CORRECT | 9 of 20 |
| mutation 2 → 1 of 13; mutation 3 → 1 of 13 | CHECKED-CORRECT | as stated, and each is the named test |
| the lib count moves **6,409 → 6,427** | CHECKED-CORRECT | 6,409 at `e99ae29b`, 6,427 at `e3e0b6e5`; the delta is exactly 13 + 5 |
| "a record with one sample and one with **a thousand**"; the message's `before.len() - 9` | CHECKED-CORRECT | the fixture builds 1,009 columns, so the message prints 1,000 |
| "with all thousand sample columns compared" | CHECKED-CORRECT | `assert_eq!(after[8..], before[8..])` compares `FORMAT` plus all 1,000 — one more than claimed, never fewer |
| "`FilterVerdict` has five values and all are words" (deviation 5) | CHECKED-CORRECT | five variants, all words |
| **"four of them tests of `write_record` written long before this step"** | **WRONG** | **seven** of the nine pre-date the commit; six call `write_record` directly, the seventh calls `write_stream`. Only `a_line_may_take_the_one_legal_tie_and_no_other` and `a_line_that_runs_backwards_is_refused_as_a_record_would_be` are new. Verified twice independently. **The claim understates its own result** |
| report §Validation, first block: "`13 passed; … 6424 filtered out`" | **WRONG** | the real line says **6429**. The report's two quoted blocks contradict each other by 5; the writer one (6,422, implying 6,442 = 6,427 + 15 ignored) is right |
| `RecordPlace` "is **exactly** the three head fields the spill entry carries" (deviation 2, commit message) | **WRONG** | two of three match; `SpillEntry` carries `position: Position`, `RecordPlace` a bare `u64` — see **M3** |
| `EMNoConv;hiddenParalog` "which the VCF grammar allows and **ng's header declares**" (`patch.rs:17`, present tense) | **WRONG at this commit** | `FILTER_DECLARATIONS` ([header.rs:398-405](../../../../src/ng/vcf/header.rs#L398-L405)) declares `PASS`, `EMNoConv`, `notPeriodic`, `tooManyAlleles`, `lowDepth`. `hiddenParalog` appears only in this step's test strings and this doc comment; spec §3.5 says C4 adds the declaration |
| the sample columns are "all but **eight** of a cohort's line" (`patch.rs:5-6`) | **WRONG** | all but **nine** |

**One mechanism claim was settled before the fan-out and is recorded here rather than filed:**
deviation 6, the `patch.rs` module doc, the commit message and a `patch/tests.rs` comment all say
ng's encoder writes an **empty** `INFO` column for a record with no annotations, "because
`info_column` joins an empty field list". It cannot —
[encode.rs:218-219](../../../../src/ng/vcf/encode.rs#L218-L219) pushes `AN=` and `DP=`
unconditionally, so the list is never empty. Confirmed at the widest input class buildable: over
all nine encoded shapes the thinnest `INFO` is `AN=0;DP=0`.

**A consequence of that worth recording as good news:** because the encoder writes neither `.` nor
`""` in `INFO`, the replace-rather-than-append branch is unreachable for any line ng emits — so the
§10 byte-identity oracle cannot be broken by it. Stripping the two keys from a *replaced* `.` would
leave an empty column where the filter-off run wrote `.`; that hazard exists in the code and not on
any reachable input.

### 7. Out of scope observations

- **[spill.rs:161-166](../../../../src/ng/run/paralog_filter/spill.rs#L161-L166) — `SpilledSample`'s
  two read-count fields document different rules, and spec §3.2 wants them to agree.** `alt_reads`
  says "`0` on every other record"; `ref_reads` says "`AD[0]`" with no qualifier. Spec §3.2 requires
  `alt_reads = 0` **and `total_reads = 0`** on every non-biallelic-SNP record, and total is
  ref + alt — so a tract that spills its real `AD[0]` gives C1 a non-zero total and a live allele
  term, which is spec §6 trap 1 arriving by a different door. Nothing fills a spill yet (pass one is
  C2), so this is latent, not a bug. It is B1's field and it bears directly on the `SpilledSample`
  shape decision waiting at Checkpoint B.
- `writer.rs:352-353` — `VcfWriteError` is public and not `#[non_exhaustive]`. Pre-existing.
- `encode.rs:262,363` — `fields.join(";")` and `out.push('\t')` are spelled as literals there too,
  so re-spelling separators is house style rather than this step's sin.

### 8. Missing tests to add now

Eight, all written and run green by the reliability agent against `e3e0b6e5`
(`ok. 29 passed` patch, `ok. 22 passed` writer). Full bodies in
[reliability_extras.md](../../../../tmp/review_2026-09-06_ng_paralog_filter_b3/reliability_extras.md).

Grouped by function under test.

**`rewrite_filter_and_info`**

1. `every_encoded_shape_including_a_three_thousand_sample_cohort_round_trips_byte_for_byte` —
   the salvaged `every_shape()` plus three shapes it misses: a record on a contig other than the
   first (every salvaged shape is `ContigId(0)`), a repeat tract whose sample is a no-call (the one
   encoded sample-column shape not exercised, where `REPCN` is `.`), and a 3,000-sample cohort
   (measured: 3,009 columns). **Catches** a splice that assumes the first contig, breaks on a
   trailing `.` in a sample column, or separates the sample columns at cohort scale.
2. `the_filter_and_info_columns_of_a_real_encoded_line_gain_exactly_the_verdict` — asserts what the
   two columns *become* on encoded lines, which neither salvaged probe does. **Catches** the two
   arguments swapped, a pre-existing verdict replaced instead of joined, the record's own `INFO`
   dropped. This is where the reachable half of trap 7 gets pinned: `notPeriodic;hiddenParalog`,
   `EMNoConv;hiddenParalog`, `lowDepth;hiddenParalog`, `tooManyAlleles;hiddenParalog`.
4. `patching_a_line_twice_adds_the_verdict_twice` — pins the absence of idempotence (Mi11).
6. `an_empty_filter_column_is_replaced_rather_than_joined_to` — **catches M5's first mutation.**
7. `a_missing_info_with_nothing_to_add_comes_back_byte_for_byte` — **catches M5's second mutation**,
   the identity law at the one `INFO` spelling the existing identity test does not use.
8. Two `proptest` properties over 10–40 columns of arbitrary non-tab bytes: nothing-to-add is the
   identity, and only two columns move (Mi14).

**`VcfWriter::write_line`**

3. `a_refused_line_leaves_the_writer_where_it_was` — **catches M6**, run against the mutant and
   failing with `left: 2, right: 1`.
5. `the_order_is_checked_against_the_place_and_not_the_line` — pins M1's boundary, so the next
   reader does not assume the bytes are checked. Replaced by a stronger test if the
   `From<&SpillEntry>` constructor becomes the only way in.

### 9. What's good

- **The splice keeps the line's bytes instead of rebuilding it from parts**
  ([patch.rs:68-79](../../../../src/ng/run/paralog_filter/patch.rs#L68-L79)) — the one decision
  that makes spec §10's byte-identity checkable rather than hoped for, and it is right.
- **The "nothing to add" case goes through the splice rather than around it** (deviation 4), so the
  byte-identity test proves the splice rather than a short circuit. Verified: the round trip holds
  over all nine encoded shapes.
- **`splitn` with `FORMAT` and the samples left as one piece**
  ([patch.rs:34-37](../../../../src/ng/run/paralog_filter/patch.rs#L34-L37)) — one split of nine at
  three thousand samples rather than three thousand of one, and a tab run inside the sample columns
  survives untouched.
- **`write_record` expressed in terms of `write_line`** — verified structurally: `sink` and `Sink`
  are both private to `writer.rs`, and the only record-bearing write is in `write_line`. No path to
  the sink skips `check_order`.
- **The tests that compare columns rather than searching for substrings**
  ([patch/tests.rs:55-69](../../../../src/ng/run/paralog_filter/patch/tests.rs#L55-L69)) — the right
  shape for a byte-identity contract, and the reason B1 is a fixture problem rather than an
  assertion problem.

### 10. Commands to re-verify

    ./scripts/dev.sh cargo test --all-features --lib "ng::run::paralog_filter::patch"
    ./scripts/dev.sh cargo test --all-features --lib "ng::vcf::writer"
    ./scripts/dev.sh cargo test --all-features --lib --bins --tests
    ./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
    ./scripts/dev.sh cargo fmt --check

### Author response convention

Address each finding by its identifier (`B1`, `M3`, `Mi7`) with one of: `fixed in <commit>` /
`disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the two open questions
in section 4 first.
