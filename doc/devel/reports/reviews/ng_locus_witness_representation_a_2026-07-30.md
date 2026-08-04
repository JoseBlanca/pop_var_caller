# Code Review: ng locus witness representation — Milestone A (the vocabulary)
**Date:** 2026-07-30
**Reviewer:** rust-code-review skill (orchestrator) + 8 category sub-agents, one isolated worktree each
**Scope:** the six rename commits `29bec44..4ad1bc1` on `ng-pileup-generator`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** a six-commit diff — 17 files, +662/−677 — implementing Milestone A of
  [the witness-representation plan](../../ng/impl_plan/locus_witness_representation.md). Six
  mechanical renames whose stated contract is that behaviour and output bytes do not change.
- **Reviewed against:** `416da4f` (code state `4ad1bc1`), branch `ng-pileup-generator`, worktree
  `pop_var_caller-ng-pileup`.
- **In-scope files:** `src/ng/locus_generation/{mod.rs, ssr.rs}`,
  `src/ng/locus_generation/pileup/{mod.rs, open_record.rs, parity.rs, tests.rs, generator.rs,
  genome_walk.rs, mock_reference.rs}`, and the seven ng example dumps.
- **Out of scope, with reason:** `src/pileup/` (production, frozen); the set types and the fold
  change (Milestones B–C, not written); `benches/psp_writer_perf.rs` (pre-existing failure).
- **Categories dispatched:** `reliability` (always), `errors` (always), `naming` (the diff *is* a
  naming change), `idiomatic` (always), `refactor_safety` (a `sed`-driven rename claiming
  byte-stability), `module_structure` (multi-file), `smells` (always), `extras` (the generators
  produce stable output, and "diff matches stated intent" applies to a six-commit series).
- **Categories skipped, with reason:** `unsafe_concurrency` — the diff touches no `unsafe`, atomic,
  lock, channel or thread; `tooling` — no `Cargo.toml`, CI or lint config changed; `defaults` — no
  default-acting value or configuration is in the diff, and the public-API rename is covered by
  `naming` and `refactor_safety`.

## 2. Verdict

**Approve-with-changes.** The milestone's central claim survives adversarial checking: three
sub-agents independently reconstructed the diff mechanically — applying the rename map to the
pre-milestone tree, stripping comments, and comparing token streams — and all three agree that **no
executable code changed semantics**. A multi-line-aware literal scanner found that exactly **four
string literals moved in the whole milestone**, of which two reach runtime output (the two dump
tools' TSV headers) and two are a compile-time assert message and a test expectation.

What the review found instead is in two groups. First, **the suite cannot see the surface the
milestone renamed** — transposing `witness_order`'s two components leaves all 275 tests green, and
transposing `ReadWitness::Partial`'s two fields at their construction site leaves all three named
anchors green. That is a Major-class gap, and it matters most for Milestones B–E, which change that
payload for real. Second, **the rename is not finished**: a public error variant still names the
removed type, ~40 prose sites still call a witness "coverage", and A6's comments-only pass left
docs and identifiers contradicting each other in the same two lines at ~40 sites.

## 3. Execution status

Each sub-agent ran the project's tooling from its **own** worktree's `scripts/dev.sh`. Results were
consistent across all eight:

| command | result |
|---|---|
| `cargo fmt --check` | clean (exit 0) |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib --bins --tests --examples --all-features` | 2,806 passed, 0 failed |
| `cargo test --release --lib ng::locus_generation` | 275 passed, 0 failed, 1 ignored |
| `cargo doc --no-deps --all-features` | 12 unresolved-link errors, 7 warnings — **all pre-existing**, none in a renamed link |

Not run: the tomato-CRAM STR oracle inside the agents' worktrees (no CRAM there); the orchestrator
ran it in the main checkout, one line differing (the STR header). `cargo test --all-targets` is red
for a pre-existing reason (`benches/psp_writer_perf.rs:386`, index out of bounds), verified on
`29bec44`.

Findings labelled "Needs verification": **0**. Every Major below was demonstrated by mutation in an
isolated worktree, with the mutation's own effect quoted.

## 4. Open questions and assumptions

1. **Is an error message "output" for the purposes of Milestone A's no-output-bytes contract?**
   Affects Mi1. `refactor_safety` says defer the rename to D4 because the `#[error]` text is
   user-visible; `errors`, `idiomatic`, `naming` and `smells` say rename now because the message
   names a type that no longer exists. *Orchestrator's resolution: rename now.* The contract's
   oracle is the generators' emitted data (the dump TSV, the locus stream); a config-rejection
   message is not that, no test asserts its text (the one that could,
   `the_constructor_rejects_a_span_no_coverage_run_could_describe`, asserts only `"65535"`), and A2
   already moved a TSV column header under the same plan step.
2. **Does "row" survive as a code word, or does the identifier rename finish now?** Affects Mi2.
   A6 renamed comments only, which the plan authorised; the consequence, which the plan did not
   anticipate, is ~40 sites where a doc and the identifier it documents disagree. *Orchestrator's
   resolution: finish it in `src/`*, leaving the two dump tools' `ObservationRow` TSV structs alone.
   The alternative the reviewers offer — say once that "row" is the emitted-table word and stop
   rewriting comments — reintroduces the two vocabularies A1–A5 removed.
3. **Should `witness_of` / `ObservationKey` / `KeyedObservation` be demoted to private?**
   Affects Mi12. Demoting compiles clean. *Orchestrator's resolution: no* — arch §2 specifies
   `pub(super) fn witness_of` and the visibility is design, not residue. Recorded, not applied.

## 5. Top 3 priorities

1. **M1/M2 — the ordering functions have no test that can fail.** Both copies of `witness_order`
   survive having their two components exchanged, with all 275 tests green. Two unit tests, both
   verified to fail under the mutation, close it.
2. **Mi10 — the census cannot see a wrong witness payload.** `classify_locus` sets `partial_witness`
   on mere presence and then skips every reconciliation, so a wrong `Partial` is its own alibi. A
   four-line per-run invariant turns the census from blind into the discriminating oracle
   Milestones B–E need.
3. **Mi1 — a public error variant still names a type that no longer exists**
   (`RecordSpanExceedsCoverageRun`), surfaced independently by five of the eight categories.

## 6. Findings

### Major

**M1: `src/ng/locus_generation/pileup/open_record.rs:263` — `witness_order` has no test that can distinguish its two components**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Exchanging the comparator's two `u16` components —
  `Partial { offset_in_locus, positions_covered } => (1, positions_covered, offset_in_locus)` —
  leaves `test result: ok. 275 passed; 0 failed`. The arm is not dead: replacing it with `panic!`
  fails four tests and prints the arm being reached with `(3,2)` and `(0,1)`. The one test named
  for the job, `parity::the_projection_orders_rows_as_the_walk_does`, sorts both sides with *this
  same function*, so it is structurally incapable of failing on any change to it.
- **Why it matters:** `witness_order` is what makes ng's emitted observation order a function of the
  observation's own identity rather than of read arrival order — a guarantee the code states at
  `open_record.rs:665-671`. An edit can silently reorder emitted evidence with nothing failing.
- **Suggested fix:** `witness_order_ranks_partials_by_offset_before_length` (§8, test 1). Verified
  by the sub-agent to pass at `416da4f` and fail under the transposition.

**M2: `src/ng/locus_generation/ssr.rs:1145` — the STR copy is untestable by the one test that names its order**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Same mutation, same green suite. `tally::tests::observed_is_sorted_by_bases_then_coverage`
  does assert left-flush before right-flush, but both its partials are built with the **same**
  `positions_covered` (`from_left(2, len 6)` and `from_right(2, len 6)`), so exchanging the
  components maps `(1,0,2)`/`(1,4,2)` to `(1,2,0)`/`(1,2,4)` and the order survives by accident.
  The fixture never varies the two components independently, which is exactly what the doc's claim
  ("offset outranks length") is about.
- **Why it matters:** On the STR path the two partial classes are different genetic constraints — a
  prefix and a suffix — and their emission order is what a cohort merge reads.
- **Suggested fix:** `tally_orders_two_partials_of_one_sequence_by_offset_before_length` (§8,
  test 2). Verified: 51 passed → 50 passed, 1 failed under the mutation.

**M3: `src/ng/locus_generation/pileup/parity.rs:1762` — the fabrication deliverable is asserted nowhere, while its twin has a floor for exactly this reason**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `measure_fabrication` computes the census headline —
  `fabricated_ref_bases += footprint.saturating_sub(positions_covered) * num_obs` — and the test
  prints it ("production credited 2,787 reads over 1,484 loci with 8,239 reference bases they never
  sequenced"). The only assertion touching the class is a *ceiling* on `fabricating_loci`
  (`parity.rs:2855`); the two numbers themselves have neither floor nor ceiling. The comment at
  `parity.rs:2864-2873` argues the stale-widen pair needs a floor precisely because "they are read
  off production's observations rather than off the classification, so a change that made every
  class-6 locus report a zero-length tail would leave the class count intact and the two numbers at
  zero". That argument applies word for word here, with no counterpart assertion.
- **Why it matters:** The fabrication triple is the generator plan's §13.2 deliverable, and it is
  derived from `positions_covered` — the field Milestones B–E replace with a set.
- **Suggested fix:** the floor in §8, test 4.

### Minor

**Mi1: `src/ng/locus_generation/pileup/generator.rs:181-192` — a public error variant, its `#[error]` text and its docs still name the removed type**
- **Categories:** errors, idiomatic, naming, refactor_safety, smells — *convergent, five of eight*
- **Confidence:** High
- **Problem:** `PileupGeneratorConfigError::RecordSpanExceedsCoverageRun` is `pub` and re-exported
  from `pileup/mod.rs:181`. Its doc (`:181`) and its rendered message (`:188`) both say "a coverage
  run", while the constant it guards (`:54`) and the function that raises it (`:117`) were both
  moved to `ReadWitness` by A2 — so the same envelope is called a `ReadWitness` run twice and a
  "coverage run" five times inside 80 lines. `grep -rn ReadCoverage src/ examples/` returns nothing,
  so this is the last `Coverage` on ng's public surface and it names a type a reader cannot find.
- **Why it matters:** It is the message an operator sees when `--max-record-span` is rejected.
- **Suggested fix:** rename to `RecordSpanExceedsWitnessRun`, change "coverage run" → "witness run"
  in the doc and the `#[error]` text, and update the four in-tree sites (`:126`, `:1288`, the
  `expect_err` at `:1287`, the re-export) plus the two test names at `:1250`/`:1363` and the prose
  at `:73`/`:606`. See open question 1 for why this lands in Milestone A rather than D4.

**Mi2: ~40 sites — A6 converted the prose but not the identifiers, so a doc and the item it documents disagree in adjacent lines**
- **Categories:** naming, idiomatic, smells — *convergent*
- **Confidence:** High
- **Problem:** The cleanest instance is a two-line field at `parity.rs:1664`:
  `/// **Class 5** — observation order: …` over `row_order: bool`. The same shape recurs at
  `open_record.rs:518` (`fn observation_rows(&self) -> Vec<KeyedObservation>`, doc "Re-derive this
  record's observations"), `parity.rs:1358` (`fn sort_rows(rows: &mut [SequenceObservation])`),
  `parity.rs:1925` (`fn row_evidence(row: &SequenceObservation, …)`, doc "The evidence on one
  observation"), `parity.rs:2244` (`type RowIdentityWithoutGroup`), `parity.rs:1719`
  (`DivergenceCensus::group_split_rows`), `tests.rs:162,198` (`reference_row` / `first_alt_row`),
  and into assertion text (`mock_reference.rs:355` comment says "an observation per read", `:365`
  message says "a row each"). The author recorded *that* four `row` identifiers survive; the
  recorded list is narrower than the reality, and the *consequence* — agreement turned into
  disagreement — is not recorded at all.
- **Why it matters:** The milestone's entire deliverable is one vocabulary. Half-converged is the
  state where a reader must hold both words and learn they mean the same thing.
- **Suggested fix:** one mechanical pass over the in-scope `src/` files (every identifier is private
  or test-only and none is printed): `observation_rows` → `keyed_observations`; `sort_rows` →
  `sort_observations`; `row_evidence` → `observation_evidence`; `rows_split_by_group` →
  `observations_split_by_group`; `RowIdentityWithoutGroup` → `ObservationIdentityWithoutGroup`;
  `DivergenceClasses::row_order` → `observation_order`; `DivergenceCensus::group_split_rows` →
  `group_split_observations`; `reference_row` → `reference_observation`; `first_alt_row` →
  `first_alt_observation`; the `_row(s)_` test names; the locals. **Leave the two dump tools'
  `ObservationRow` structs and their `rows` fields alone** — those are TSV rows of a real table.

**Mi3: `src/ng/locus_generation/pileup/open_record.rs:110-137` — `RecordWitness`'s doc states two things that are now false, and its field name spells the variant A4 deleted**
- **Categories:** naming, idiomatic, smells — *convergent*
- **Confidence:** High
- **Problem:** The doc says "only until B2. `finalise` still returns production's [`PileupRecord`],
  which has nowhere to put a [`ReadWitness`]". `finalise`'s signature at `:593` is
  `-> (SampleLocusObservations, RecordWitness)`; B2 landed. Separately the field
  `reads_partially_observed` (`:127`) counts the `ReadWitness::Partial` arm, so its name says
  `Observed` where the code says `Partial`. And the local `witness: RecordWitness` in `finalise`
  sits one line from `witness_of(...) -> ReadWitness`, so the word now names three things in the
  file (the per-read classification, the per-record tally, and the `witnessed` extent).
- **Why it matters:** A reader is told `finalise` returns production's type — the most load-bearing
  fact about the walk's current shape — and it is wrong.
- **Suggested fix:** rewrite the doc; rename the field to `reads_partial`; rename the type to
  `RecordWitnessCounts` and the local to `witness_counts`. **The type rename is already sanctioned**
  — the [arch doc](../../ng/arch/locus_witness_representation.md) §3's boy-scout note asks for
  exactly it, "when the file is next touched".

**Mi4: ~40 prose sites across the in-scope files — "coverage" still names the concept now called a witness, in a codebase where "coverage" also means depth**
- **Categories:** naming, idiomatic, smells, errors, extras, reliability — *convergent, six of eight*
- **Confidence:** High
- **Problem:** Exemplars: `mod.rs:140` "must aggregate over **coverage** *and* group" (the axis is
  `read_witness`); `mod.rs:970` "at a given coverage"; `open_record.rs:228` "coverage and group are
  facts about a read"; `open_record.rs:662` the `debug_assert` message "every folded read resolves
  to exactly one **coverage** class", eight lines below a match on `ReadWitness`;
  `parity.rs:1966-1971` four times in one doc. The genuine depth sense must stay: `mod.rs:46`,
  `chain_id_allocator.rs:48-59`, `parity.rs`'s `windowed_coverage`, `ssr.rs`'s "zero coverage".
- **Why it matters:** A2's own rationale is that "coverage" misreads as depth to a geneticist. The
  prose that does the misreading is what stayed, and the type name that used to disambiguate it is
  gone.
- **Suggested fix:** one sweep, witness sense → "witness", depth sense untouched.

**Mi5: `src/ng/locus_generation/mod.rs:130-143` — `SequenceObservation`'s doc now contradicts itself**
- **Categories:** naming, smells
- **Confidence:** High
- **Problem:** "One distinct sequence the reads showed at a locus" … then "**This is a table of
  observations, not a table of sequences.**" Pre-A6 the contrast was "a table of **cells**, not a
  table of sequences" and it worked, because "cell" was the corrective noun. A6 kept the rhetorical
  shape and removed the contrast — the type is now literally `SequenceObservation`.
- **Why it matters:** This paragraph is the public warning against treating an entry as an allele,
  and it now argues against its own first sentence.
- **Suggested fix:** replace the contrast with the thing it was trying to say — one entry is not one
  allele; the identity has three axes.

**Mi6: `examples/ng_generic_loci_dump.rs:267` — the generic dump's TSV header is pinned by no test, and the milestone moved it**
- **Categories:** refactor_safety, extras — *convergent*
- **Confidence:** High
- **Problem:** Two runtime output literals moved, not one: the STR dump's header *and* the generic
  dump's. The plan authorises both ("the `read_coverage` column in **both** dump tools' TSV
  output"), but the STR dump is the only one with a pin
  (`render_emits_the_spec_9_header_and_tsv_rows`) and the only one the tomato-CRAM oracle covers.
  The generic dump's `render()` is exercised only by two self-comparisons (`first == second`,
  `huge == whole`), which pass whatever the header says. A sub-agent confirmed the change was real
  only by building both commits and diffing the rendered output itself.
- **Why it matters:** Milestones B–E rewrite exactly this rendering; today a dropped column or a
  reordered header there passes CI silently.
- **Suggested fix:** mirror the STR dump's header assertion (§8, test 5).

**Mi7: `src/ng/locus_generation/pileup/tests.rs:392` — the substitution produced a sentence that says nothing**
- **Categories:** smells
- **Confidence:** High
- **Problem:** `"one read, one observation in its observation"`. Pre-milestone: `"one read, one
  observation in its row"`, where "observation" meant a read's contribution and "row" the table
  entry. Six lines above, `tests.rs:386` still says "Rows are derived from `folded_reads`", so one
  comment block carries both vocabularies and the retired one is the only sentence in it that
  parses.
- **Suggested fix:** "one read contributes `num_obs: 1` to the observation it lands in".

**Mi8: `src/ng/locus_generation/ssr.rs:1215` — A6's own headline example did not move, and it is the last "cell" in the tree**
- **Categories:** smells
- **Confidence:** High — *orchestrator re-verified: the line is byte-identical to `29bec44`*
- **Problem:** The A6 commit message opens by naming `ssr.rs`'s "one row per (allele, read group)
  cell" as a nickname that "moves now". It did not: it is an assertion message, and A6 restricted
  itself to comments. It is the only "cell" in the observation sense left in
  `src/ng/locus_generation/`.
- **Why it matters:** A commit message claiming a specific string moved when it did not is a false
  record for whoever audits the milestone's completeness.
- **Suggested fix:** change the message; correct the record in the fix report.

**Mi9: `src/ng/locus_generation/pileup/parity.rs:2143,2162` — "observation" now collides with `num_obs`, and the assertion reports a read count as an observation count**
- **Categories:** smells
- **Confidence:** High
- **Problem:** `ours_obs` is `observations.iter().map(|row| row.num_obs).sum()` — a **read** count —
  and the message reports it as "ng emitted {ours_obs} observations". `observations.len()` is the
  observation count; on any locus where two reads share an identity the two differ. Both lines are
  unchanged from `29bec44` and were unambiguous then, when entries were called rows.
- **Suggested fix:** "ng's observations carry {ours_obs} reads"; and in the doc above, "contributing
  exactly one to that bucket's `num_obs`".

**Mi10: `src/ng/locus_generation/pileup/parity.rs:2067-2129` — none of the three named anchors discriminates on the surface this milestone renamed**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Transposing `ReadWitness::Partial`'s two fields at `witness_of`'s construction site
  fails five unit tests and leaves all three anchors green. Three mechanisms, each verifiable:
  the anchor's fixture yields no partial witness at all (216,203 of 216,203 loci complete, per its
  own doc); the determinism digest compares two children of the same binary, so deterministic
  corruption hashes identically; and `classify_locus` sets `partial_witness` on mere *presence*,
  after which it skips the per-`bases` reconciliation, skips the chain-id equality, and accepts any
  difference via `exact || classes.any()` — **a wrong `Partial` payload is its own alibi.**
- **Why it matters:** The milestone's warrant is that the existing suite proves the rename changed
  nothing. It does — but the proof comes from `open_record.rs`'s and `mod.rs`'s unit tests, not from
  the anchors the plan cites, and the anchors will be equally silent when B–E change the payload for
  real.
- **Suggested fix:** the per-run invariant in §8, test 3 — verified to fail the census under the
  transposition with a named locus, and green unmutated.

**Mi11: `src/ng/locus_generation/pileup/open_record.rs:263` + `src/ng/locus_generation/ssr.rs:1145` — two `witness_order` functions with byte-identical bodies, in a tree that argues in prose against exactly that**
- **Categories:** module_structure, smells — *convergent*
- **Confidence:** High
- **Problem:** The duplication pre-dates the milestone (both were `coverage_order`), but A3 renamed
  both and the collision is now visible by name. `open_record.rs:257-259` justifies withholding an
  `Ord` impl because it "would export **this file's** sorting convention to every other consumer" —
  a claim the `ssr.rs` copy refutes: two independent files want the identical convention, so it is
  the type's. A sub-agent verified that deriving `PartialOrd, Ord` on `ReadWitness` and replacing
  both call sites gives 275 passed, 0 failed, parity tests included.
- **Why it matters:** B replaces `Partial`'s payload; both copies must be rewritten, and the
  ordering semantics each is rewritten to is a free choice per site — the two paths' emission orders
  can diverge silently while both compile.
- **Resolution: deferred to Milestone B**, where `witness.rs` is created and both types move. B1's
  move should absorb it rather than carry two copies across. Recorded in the plan.

**Mi12: `src/ng/locus_generation/pileup/open_record.rs:186,230,239` — three renamed items are `pub(super)` with no consumer outside their own file**
- **Categories:** module_structure
- **Confidence:** High
- **Problem:** `witness_of`, `ObservationKey` and `KeyedObservation` are `pub(super)`; demoting all
  three to private leaves `cargo clippy --all-targets --all-features -- -D warnings` clean. The
  neighbouring `witness_order` is the contrast — its `pub(super)` is documented and load-bearing for
  `parity.rs`.
- **Resolution: disputed.** Arch §2 specifies `pub(super) fn witness_of` in the interface this plan
  is implementing, and C2 changes its signature there. The visibility is a design decision, not
  rename residue. Recorded so it is not re-found.

**Mi13: `src/ng/locus_generation/pileup/{open_record.rs:27, genome_walk.rs:28, parity.rs:122}` — the witness vocabulary is imported by `super::super::` in three files and crate-absolute in two**
- **Categories:** module_structure
- **Confidence:** High
- **Problem:** A1/A2 rewrote exactly these import lines. In `open_record.rs` the same spelling
  resolves to two different modules: `locus_generation` at file top, `pileup` inside the nested test
  module. `generator.rs:19` and `tests.rs:33` already use the crate-absolute path for the same types.
- **Resolution: deferred to Milestone B**, which moves the types into `witness.rs`; converting the
  three import sites then makes `grep crate::ng::locus_generation::witness` answer "who depends on
  the witness vocabulary" in one command.

**Mi14: the commit messages' accounting is wrong in four places**
- **Categories:** extras, refactor_safety, smells — *convergent*
- **Confidence:** High
- **Problem:** (a) A2 says "No expectation was edited"; it hand-edited the header literal in
  `ng_ssr_loci_dump`'s `render_emits_the_spec_9_header_and_tsv_rows` — a string, not an identifier,
  so no rename tool moved it. This matters because the plan's tripwire is literally "any test that
  changes expectations did more than rename", and the sentence disarms the tripwire instead of
  naming the one benign instance. (b) The implementation report says "The one output change in the
  whole milestone is the STR dump's column header"; two headers moved. (c) A6 says "five sentences
  were rewritten" and names three, and says "two article fixes" where a mechanical count gives ~23.
  (d) A6 claims `ssr.rs:1215` moved (see Mi8).
- **Suggested fix:** correct all four in the implementation report and the fix report. Do not
  rewrite the commits.

**Mi15: `examples/{ng_ssr_loci_dump.rs:163, ng_ssr_cohort_stutter.rs:154, ng_ssr_aligner_bakeoff.rs:198}` — three copies of `witness_label`, identical down to a shared seven-line comment, already drifted**
- **Categories:** smells
- **Confidence:** High
- **Problem:** Same signature, same match shape, same comment verbatim; they differ only in label
  strings, and that difference is already incoherent — one emits `partial:left`/`partial:right`,
  the other two emit `partial_left`/`partial_right` but keep `partial:interior`. Pre-existing
  (verified at `29bec44`); A2 renamed all three, which is what made it visible.
- **Resolution: deferred to D4**, which rewrites the labels when the dumps start printing the set.
  Recorded in the plan so the `partial:interior` inconsistency is decided rather than inherited.

### Nits

Grouped, per the volume guidance — all mechanical:

- **Article agreement.** The A4 substitution produced "an `Partial`" at `mod.rs:61`, `mod.rs:149`,
  `parity.rs:1646` and `examples/ng_generic_loci_dump.rs:1151`. Two are public rustdoc.
- **Comment wrap.** "observation" is eight characters longer than "row" and nothing re-flowed;
  comment lines over 96 columns went **46 → 99** across the five main files (`cargo fmt` does not
  reflow comments). Best done in the same pass as Mi4 so the churn lands once.
- **Stale locals and helpers.** `let coverage = witness_of(…)` at `open_record.rs:2947,2963,2981`
  (inside the very tests A3 renamed); the fixture `same_bases_different_coverage`; `let coverages`
  at `ng_generic_loci_dump.rs:1171`; `ssr.rs:651`'s field `witness` where every other site spells
  the same concept `read_witness`.
- **Test names spelling retired words.** `a_complete_witness_becomes_observed_when_the_record_widens_under_it`
  names the variant A4 deleted; `observed_is_sorted_by_bases_then_coverage`;
  `a_span_one_past_the_ceiling_is_where_a_coverage_run_starts_lying`;
  `the_constructor_rejects_a_span_no_coverage_run_could_describe`.
- **~25 assertion messages still say "row"/"cell"** — the message-string half of A6's
  comments-only limit.
- **The per-step gate excludes `--benches`**, so a rename breaking a bench's *compile* would not be
  caught. Checked: no bench references any renamed item, so there is no live breakage.

## 7. Out of scope observations

- **`doc/devel/ng/arch/locus_generation.md`, `arch/locus_generation_pileup.md`,
  `arch/locus_generation_ssr.md`, `arch/module_layout.md`, `arch/ng_step_interfaces.md`,
  `doc/devel/ng/README.md`** carry ~19 combined occurrences of `ObservedSequence` / `ReadCoverage` /
  `coverage_of` / `observed_sequences`. The implementation report records only the first. Follow-up:
  a docs pass, not part of this plan.
- **`generator.rs` and `pileup/tests.rs` use `..Config::default()` in 14 test-fixture struct
  literals** — the sites that will silently absorb a new config field. Pre-existing, untouched by
  this diff; belongs to a whole-file review of the generator.
- **`mod.rs:1638` cites `ssr::tally::tests::an_expanded_allele_merges_the_two_sides`**, but the test
  is `…_into_one_row`. Pre-existing in both spellings at `29bec44`.
- **`mod.rs:1658` asserts a property of `witness_label`**, an identifier that exists only in three
  `examples/` binaries and nowhere in the crate.
- **`benches/ng_generic_pileup_perf.rs` carries no regression threshold.** Its own header records
  that performance is parked on the owner's instruction until correctness is done.

## 8. Missing tests to add now

1. **`witness_order_ranks_partials_by_offset_before_length`** — `open_record.rs`, `mod tests`.
   Two `Partial` runs whose components vary in opposite directions (`{0,9}` vs `{4,2}`), never
   produced by any existing fixture. Catches: the two components exchanged, or the
   `Complete`/`Partial` ranks exchanged. *Verified: passes at `416da4f`, fails under the mutation.*
2. **`tally_orders_two_partials_of_one_sequence_by_offset_before_length`** — `ssr.rs`,
   `mod tally::tests`. One sequence witnessed by a **long** left-flush run and a **short**
   right-flush run. Catches the STR copy's components exchanged. *Verified: 51 → 50 passed, 1
   failed under the mutation.*
3. **The census's per-run invariant** — `parity.rs`, inside `classify_locus`, beside the
   region/reference-bases assertion so it runs at every locus of both census passes: every `Partial`
   has at least one position and lies inside its own locus. Catches any wrong payload from
   `witness_of`, including the transposition all three anchors sleep through. *Verified: green
   unmutated; under the transposition the census fails naming the locus.* It is an **assertion, not
   a divergence class** — a run outside its own locus is ng being wrong on its own terms, not a
   difference from production.
4. **`the_fabrication_deliverable_is_positive_wherever_the_class_fires`** — `parity.rs`, beside the
   stale-widen floor: where class 1 fires, `fabricated_reads > 0` and
   `fabricated_ref_bases >= fabricated_reads`. Closes the asymmetry the code itself argues for.
   Confidence Medium on its power — it is a floor, so it does not catch the transposition; test 3
   is what does.
5. **The generic dump's header assertion** — `examples/ng_generic_loci_dump.rs`, mirroring
   `ng_ssr_loci_dump`'s `render_emits_the_spec_9_header_and_tsv_rows`.

Deliberately **not** added: a test for `witness_of` on an extent lying wholly past the footprint.
The function opens with a `debug_assert!` that the input is unreachable, so an assertion on the
returned value passes in release and fails in debug, and `#[should_panic]` does the reverse — either
spelling is profile-dependent. The honest resolutions are to gate two tests explicitly by profile,
or to decide the disjoint case is unreachable and delete the clamp. **Author decision, recorded as
an open item.**

## 9. What's good

- **The renames themselves are the right ones and are justified in place.** `SequenceObservation`
  matches production's own `AlleleObservation` shape (`src/pileup_record.rs:138`) where
  `ObservedSequence` did not, and `ReadWitness::Partial` carries a four-line record of why it is
  not `Observed` (`mod.rs:216-222`).
- **The two blanket substitutions that would have been silent damage were caught and handled**:
  `ssr.rs`'s unrelated `Classified::Observed` (9 sites intact) and `generator.rs`'s
  `std::cell::Cell` comments (reverted whole, then hand-edited). Every English idiom survived —
  "two `u16`s in a row", "twice in a row", "the debug row".
- **The step-per-commit discipline made the mechanical verification possible.** Three sub-agents
  independently reconstructed the diff by applying the rename map to `29bec44`; that only works
  because no behaviour change is mixed in.
- **`open_record.rs:665-671`'s comment records that the row sort was mutation-tested both ways** —
  deleting it leaves the determinism test green, deleting `ids.sort_unstable()` fails it. That is
  the kind of note the reliability findings above wish existed for `witness_order`.

## 10. Commands to re-verify

```
# from THIS worktree's wrapper — the main repo's copy builds the main worktree
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo fmt --check
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo test --lib --bins --tests --examples --all-features
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo test --release --lib ng::locus_generation
```

New with this review — each must fail under its own mutation, not merely pass:

```
… cargo test --release --lib ng::locus_generation::pileup::open_record::tests::witness_order_ranks_partials_by_offset_before_length
… cargo test --release --lib ng::locus_generation::ssr::tally::tests::tally_orders_two_partials_of_one_sequence_by_offset_before_length
… cargo test --release --lib ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes
… cargo test --release --example ng_generic_loci_dump
```

The external oracle, unchanged:

```
DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks/tomato1 \
  … cargo run --release --example ng_ssr_loci_dump -- \
  /Users/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
  /Users/jose/devel/pop_var_caller/benchmarks/tomato1/crams/SRR7279503.p1.bench.cram SL4.0ch01
```
