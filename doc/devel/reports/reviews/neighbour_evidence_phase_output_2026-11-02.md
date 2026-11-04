# Code review: writing the phase the neighbour-evidence pass estimates (2026-11-02)

*Merged report of a review by the `rust-code-review` skill: five sub-agents, each in its own
worktree, over nine checklists and the check of the diff's own numbers. Per-category files are in
the `neighbour-pass` worktree's `tmp/review_2026-11-02_phase-output/`. Each finding below carries
its disposition; the fixes are commit `584f4d48` and the document corrections after it.*

## 1. Scope

- **Reviewed:** the diff `d31e4477..19102f1e` on branch `neighbour-pass`, the build of
  [`neighbour_evidence_phase_output.md`](../../ng/spec/neighbour_evidence_phase_output.md).
- **In scope:** `src/neighbour_evidence/phase.rs` (new) and the diffs of `refine.rs`, `model.rs`,
  `random.rs`, `stage.rs`, `line.rs`, `test_cohort.rs` and their tests, `src/vcf/header.rs`,
  `src/cli/{refine_genotypes,calling_run,call_from_alignments,call_from_psps}.rs`; the commit
  messages, the implementation report and the user page for the numbers check.
- **Out of scope:** the rest of the crate.
- **Categories:** reliability (mutation testing); float portability; errors, defaults and
  refactor safety; naming, idiomatic code and module structure; smells, extras and the numbers
  check. `unsafe_concurrency` was not dispatched: the diff adds no lock, atomic or `unsafe`, only
  work inside the existing per-sample `rayon` closures. `tooling` was not dispatched: no
  `Cargo.toml` change.

## 2. Verdict

**Approve-with-changes**, the changes made. No finding was about a wrong genotype or a wrong
phase in the output; the one Blocker was a test that could not see a wrong phase.

## 3. Execution status

Run natively on macOS arm64. `cargo fmt --check`: differences in three files this diff does not
touch. `cargo clippy --lib --tests`: clean (`--all-targets --all-features` fails on
`tests/ng_calling_loop_allocation.rs`, untouched and failing before). `cargo test --lib`: 5,086
passed before the fixes, 5,102 after. `cargo doc --no-deps`: clean. Reliability: 33 mutations
run, 13 killed, 10 survived with a proven change of behaviour, 10 changed no behaviour.

## 4. Open questions and assumptions

- The written phase's bits were recorded on macOS arm64 only (`0xd508e07338527b7f` over 2,574
  phased columns); Linux x86-64 must confirm them, as for the posteriors' hash.

## 5. Top 3 priorities

1. B1 — the join test could not see a reversed join. **Fixed.**
2. M1 — more than 8 rounds panicked in the middle of a window. **Fixed.**
3. M4 — nothing pinned the written phase's bits. **Fixed.**

## 6. Findings

### Blocker

- **B1. `blocks_join_across_windows` could not catch a join written reversed** (reliability).
  Comparing copy 1's allele instead of copy 0's at the join flipped all 29 links that cross a
  window boundary, and the test still passed: it counted the crossing flips and asserted only the
  total under 1%. **Fixed:** it asserts at least 10 crossing links and at most a tenth of them
  wrong (when written, 29 and 0).

### Major

- **M1. `Phasing::new` accepted round counts whose paths overflow the 64-bit link masks**
  (reliability, errors, idiomatic, naming, smells — convergent). 9 rounds draw 72 paths and
  panicked at `record_path`, inside a window. No command-line value reached it. **Fixed:**
  `Phasing::new` refuses it with the reason, `MAX_PATHS` names the limit, and a `const` assertion
  ties the defaults to it; test `more_rounds_than_the_link_test_holds_are_refused`.
- **M2. The "leans at both records" half of the link test was untested** (reliability): removing
  it let 39 wrong links in 7,466 through on the stage fixture, against 9 in 5,494. **Fixed:**
  `a_link_breaks_where_one_path_does_not_lean`.
- **M3. The join rule's choice of records was not pinned** (reliability): taking the first three
  shared records, or the window's first phased record rather than its first core one, passed.
  **Fixed:** `a_join_reads_the_records_nearest_the_core`,
  `a_join_reads_the_block_of_the_first_core_record`.
- **M4. No test pinned the written phase, and none the extra paths' stream** (reliability, float
  portability): drawing them from the round's own stream changed 2,574 phased columns to 2,596 and
  passed. **Fixed:** `a_phased_windows_phase_is_the_same_bits_everywhere` hashes every written
  phase; unchanged by the fixes.
- **M5. Four changes to the alignment survived** (reliability). **Fixed in part:** a tie test and
  a test of a round that does not lean kill two; the phase digest (M4) covers a changed
  `SWITCH_COST`. The back-pointer tie rule was not re-tested.
- **M6. The cross-window state was exercised only where it made no difference** (reliability):
  four mutations of `remember_last_blocks` changed nothing, and the older one-against-four-threads
  test writes no phase. **Fixed in part:** `the_phase_is_the_same_at_any_thread_count` (1, 4, 8
  threads, a run that phases) and `a_phase_set_carries_through_short_windows` (cores of 20–36 kb
  under 60-kb margins; 95 crossing links, none wrong; a phase set reaching into 3 cores). **Still
  open:** not extending the last block when it continues was re-tried against the new test and
  changed nothing (5,693 links, 59 wrong, 95 crossings either way) — it matters only when a core
  holds fewer than 3 of the block's records, which no fixture builds.
- **M7. An `expect` without its `PANIC-FREE` note in `written_in_core`** (errors). The invariant
  held. **Fixed:** the join carries each member's phase, so no `expect` is left.

### Minor

| | finding | categories | disposition |
|---|---|---|---|
| Mi1 | the path-agreement share a bare `f64` checked in three places, one of them bypassable through a `pub` field | errors, idiomatic, smells, module structure | **fixed**: `PathAgreement`, checked once |
| Mi2 | `first`/`second` meant three things; `block` held the `PS`; functions named as values | naming | **fixed**: `copy_0_allele`, `copy_1_allele`, `phase_set`, `phase_heterozygotes`, `join_and_name_core_blocks`, `align`, `SHARED_RECORDS_TO_JOIN` |
| Mi3 | the last written block a 4-number tuple, in three shapes across stage, refine, phase | idiomatic, smells, module structure | **fixed in part**: a named struct; two shapes left (the stage's run-wide sites, the window's indices) |
| Mi4 | `consensus()` allocated two vectors a site a round | idiomatic, extras | **fixed** |
| Mi5 | `record_path` recomputed each pair's index | extras | **fixed**: a running index |
| Mi6 | GT/PS declared by two rules that disagreed on where PS goes | module structure | **fixed**: one function, PS beside GT |
| Mi7 | the new GT constant took the doc comments of `FORMAT_DECLARATIONS` and of `GP` | naming, smells, reliability, errors | **fixed** |
| Mi8 | `resize` cut a column longer than its `FORMAT` | errors | **fixed** |
| Mi9 | a last block's record checked "not in the core" rather than "in the left margin" | errors | **fixed** |
| Mi10 | `AlreadyPhased` gave a record's first 80 characters; the refusal comes after the output is opened | errors | **fixed in part**: names `contig:position`; still refused while streaming |
| Mi11 | the header recorded neither the join count nor the alignment costs | defaults | **fixed**: `join_records`, `alignment_costs` |
| Mi12 | defaults filled in without saying so in `Phasing::new`, `rounds_for`, `decide` | defaults | **fixed in part**: `Phasing::new`'s doc and the decision's field; `rounds_for` refers to `new` |
| Mi13 | `sites` and `phases` zipped, so a length mismatch would drop outcomes silently | refactor safety | **fixed**: asserted |
| Mi14 | `same[0]` would panic with a join count of 0 | refactor safety | **fixed**: a slice match |
| Mi15 | `ceil` → `floor` in the paths needed survived | reliability | **fixed**: `a_share_asks_for_whole_paths_rounded_up` |
| Mi16 | the option's parser and the share's bounds untested | reliability, errors | **fixed**: two tests |
| Mi17 | no test of a phase on an uncalled column (spec §7.4) | reliability, extras | **fixed**: `a_phase_on_an_uncalled_column_is_not_written` |
| Mi18 | the window test's 1% ceiling sits far above its measured 0 in 3,341 | reliability | open: left; the digest (M4) now pins the exact output |
| Mi19 | `final_solution` about 150 lines, `phased` 77 with `previous` tested twice | smells | **fixed in part**: the link loop rewritten; `final_solution` left |
| Mi20 | spec §2.3 said the extra draws share one pass a record; the code draws each separately | extras, claims | **fixed** in the spec |
| Mi21 | nothing pinned the written phase across platforms | float portability | **fixed in part**: the digest (M4); Linux x86-64 run open |
| Mi22 | new `pub` items with no caller outside the module | idiomatic | open, left as is |

### Wrong numbers (the numbers check)

Each corrected where it stands in a document or a doc comment; the commit messages of `c2b2c13e`
and `19102f1e` keep the old figures.

| claim | was | is |
|---|---|---|
| joins made on four cohorts | 57 | 53 |
| joins added with 2 shared records | 13 | 11 |
| "Rust phases more on every cohort" | every | seven of nine (fewer at 10 outbred samples, none at 10 inbred) |
| cohorts where Rust's links are wrong more often | three, all under 1% | four, one far over (10 outbred: 7 of 13) |
| the line test's example | `0|2` written `0|1` | `2|0` written `1|0`, repeat counts `9,7` |
| links all paths agree on, wrong | 0.05% on one cohort | 0.06%; and 37% at 10 outbred samples |
| spec §2.3, how the extra draws are made | one shared pass | a backward draw each |

## 7. Out of scope observations

- `tests/ng_calling_loop_allocation.rs` does not compile with `--all-features`.

## 8. Missing tests to add now

The tests listed under the findings were added. Still missing: a fixture where a core holds fewer
than three of a continuing block's records (M6), and the back-pointer tie in the alignment (M5).

## 9. What's good

- The dependency direction is clean: `vcf/header.rs` imports nothing from the pass, `phase.rs`
  only `model::HaplotypePair`.
- Every float decision is portable: no transcendental call in the new code, per-sample sums in a
  fixed order, ties broken by explicit rules, the paths-needed count exact for any share.
- The extra paths change neither the genotypes nor anything a round passes on, and a test proves
  it bit for bit.

## 10. Commands to re-verify

```
cargo test --lib neighbour_evidence
cargo test --lib -- cli::calling_run cli::refine_genotypes vcf::
cargo clippy --lib --tests
```
