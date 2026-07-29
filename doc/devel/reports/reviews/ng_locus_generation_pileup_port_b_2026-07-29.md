# Code Review: ng_locus_generation_pileup_port_b
**Date:** 2026-07-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Milestone B of the ng generic locus generator port — the stage-1 parity oracle
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `8a7368a..fcf0c3e` on branch `ng-pileup-generator` — three commits,
  `0a77b2e` (B1), `f39f125` (B2), `fcf0c3e` (B3). 1,044 insertions, almost all of them one new
  `#[cfg(test)]` file.
- **In-scope files:**
  - [pileup/parity.rs](../../../../src/ng/locus_generation/pileup/parity.rs) (new — the harness)
  - [prepared_read.rs](../../../../src/ng/read/prepared_read.rs) (changed — `into_production`)
  - [pileup/mod.rs](../../../../src/ng/locus_generation/pileup/mod.rs) (changed — one `mod parity;`)
- **Deliberately out of scope:** the contents of the eight copied walker files (a verbatim copy of
  frozen production code, enforced by `copy_fidelity.rs`) — though two reviewers mutated them
  *temporarily* and restored them, verified byte-exact by `copy_fidelity` passing afterwards. Also
  `src/pileup/`, `src/psp/`, `src/var_calling/`, `src/vcf/` (frozen) and the pre-existing bench
  failure.
- **Categories dispatched (6):** `reliability` (the highest-stakes category here — this is the
  baseline plan 3 measures against, and it cannot be reconstructed later), a dedicated
  **generator-correctness** pass (the differential is worth exactly what its generator reaches),
  `errors`, `naming`, `idiomatic`, `smells`.

**Two reviewers went beyond reading.** The reliability pass **re-ran two of B2's five mutations** and
reproduced the recorded table to the seed, the case index and the stream item — then ran a *sixth*
mutation the table did not cover and found a hole. The generator pass **re-implemented the generator
in a scratchpad**, verified fidelity against the harness's own printed counter (4542 = 4542), and
measured the population over 18,077 generated reads. Neither was asked for; both are why this review
found what it did.

## 2. Verdict

**Approve-with-changes.** The oracle is sound in its core claim and was demonstrably able to fail —
but two of its stated claims were stronger than its code, and the generator had a real coverage hole.
All applied.

## 3. Execution status

Verbatim output in `tmp/review_2026-07-29_ng-pileup-port-b/verification.txt`.

| command | exit | result |
|---|---|---|
| `cargo fmt --all -- --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | no diagnostics |
| `cargo test --all-targets --all-features` | 101 | lib suite **2644 passed / 0 failed / 5 ignored** at review time; the sole failure is the pre-existing `benches/psp_writer_perf.rs:386` panic |

`cargo doc --no-deps` and `cargo audit` not run (known-red / not installed). Findings labelled "Needs
verification": **0** — every number below was measured.

## 4. Open questions and assumptions

1. **The production defect the harness found** (`apply_events_to_ref_into`'s reachable
   `debug_assert!`, 4.2% of generated cases) is recorded, not fixed, because production is frozen.
   Does it want a research note of its own? Carried to Checkpoint B.

## 5. Top 3 priorities

1. **M1** — the panic *cause* was discarded, so ng panicking for an unrelated reason passed. Proven
   by mutation, with the whole suite green.
2. **M2** — the generator emitted only 5 of 9 `CigarOp` variants; `=`/`X`/`H`/`P` never appeared.
3. **M3** — the real-data differential had no floor on records compared, so a walk dying on read one
   printed "1 records compared, zero divergences" and passed.

## 6. Findings

### Blocker

None.

### Major

**M1: parity.rs — `WalkOutcome.panicked` was a `bool`, so two walkers panicking for unrelated reasons were recorded as agreeing**
**Categories:** reliability, errors (convergent) · **Confidence:** High · **Applied**

`catch_unwind(..).is_err()` discarded the payload. The reviewer replaced ng's copy of production's
reachable `debug_assert!` with `panic!("NG PANICS FOR A COMPLETELY DIFFERENT REASON…")` and **the
entire parity module stayed green** — including `both_walkers_panic_on_a_deletion_anchored_before_its_record`,
whose own doc claimed to check that ng "reach[es] the same precondition on the same input".
`open_record.rs` alone carries eight distinct `debug_assert!`s. The mutation is not contrived: it
converts a debug-only assert into a release-live panic, the transcription slip that costs most.

**Fix applied:** `WalkOutcome` carries `panic_message: Option<String>`, downcast from the payload;
`assert_same_walk` compares the messages and does so **first**, so a panic divergence is diagnosed as
one rather than as a record-count mismatch. `both_walkers_panic_on_…` now asserts the message
contains `apply_events_to_ref_into: event anchor`, pinning the defect **by cause**. Re-verified: the
mutation now fails at seed 0 case 15 with "the two walkers did not stop the same way".

**M2: parity.rs — the generator emitted only 5 of the 9 `CigarOp` variants**
**Categories:** generator · **Confidence:** High · **Applied**

`SeqMatch` (`=`), `SeqMismatch` (`X`), `HardClip` and `Padding` never appeared — 0 of 18,077 measured
reads. `=`/`X` are what minimap2 `--eqx` and DRAGEN emit and share the `Match` arm at four cursor
sites; `H`/`P` consume neither axis and are what probe the offset-table walk. A transcription slip
dropping `SeqMismatch` from one arm would have stayed green over any soak.

**Fix applied:** the match block now draws `M`/`=`/`X`; a leading `HardClip` on 1-in-6 reads and an
interior `Padding` were added. Coverage re-measured; the differential still passes.

**M3: parity.rs — the real-data differential had no floor on records compared, and no check that any record was `Ok`**
**Categories:** reliability · **Confidence:** High · **Applied**

The only floor was `!ng_reads.is_empty()` — one read. Every `WalkerError` is fatal and terminal, so a
walk dying on its first read yields one identical `Err` on each side, `assert_same_walk` agrees, and
the test prints `"1 records compared, zero divergences"` and passes. The synthetic sibling has exactly
this floor; the real-data one did not — and it is the only evidence in the milestone that the two
walkers agree on real data, hand-run, with the numbers recorded in a doc comment.

**Fix applied:** asserts no `Err` in the stream, and `ok_records * 4 > prepared_reads`. The message
now reports records *and* prepared reads, so a thin run is visible in its own output.

**M4: parity.rs — `reads_with_live_adaptor_boundary` was exactly `adaptor_boundary.is_some()`**
**Categories:** reliability, generator (convergent) · **Confidence:** High · **Applied**

Both arms returned `true`, so the "falls inside their own span" claim was never checked. Measured:
`counted_live = 4542`, `boundary_set = 4542` (identical), and of those **126 (2.8%) silence nothing**
— `base_in_adaptor` is consulted only at Match-emit sites, so a forward boundary landing in a
trailing `D`/`N` tail is inert. So `assert!(adaptor_boundaries > 0)` could not fail while
`one_in(4)` fired once in 1600 cases: **the one assertion written to prevent "a test that cannot
fail" was one.**

**Fix applied:** liveness is computed the way the cursor computes it, over the positions a
`Match`/`=`/`X` will actually emit, direction-aware. And the counter's real claim is now checked end
to end by a new test, `the_adaptor_filter_changes_the_records_the_walk_emits`: clearing every
boundary must change some case's records.

**M5: parity.rs — the differential never compared a single `Err` item**
**Categories:** reliability · **Confidence:** High · **Applied**

Measured `0 walker errors` over 1600 cases. The generator is built to stay in bounds and in order, so
`assert_eq!` only ever saw `Ok` and the `map_err` machinery comparing two nominally distinct
`WalkerError` types was dead — while spec §3 states the claim as the two
`Result<PileupRecord, WalkerError>` streams being equal element for element.

**Fix applied:** `both_walkers_report_the_same_error_on_the_same_malformed_input`, four fixtures
(out-of-order, zero ref span, CIGAR/`seq` mismatch, `seq`/`bq` mismatch), each **required to reach
its error** — which immediately earned its keep: the zero-ref-span fixture did not, because the check
is `alignment_end < alignment_start` rather than "consumes no reference". Separately, the contig-end
placement from Mi3 now produces 134 genuine walker errors in the main differential too.

### Minor

**Mi1 — the fixup for a zero-reference-span read was dead and latently unsound.** Unreachable
(needs `ref_pos > 148`, max 57), and it clamped the base *index* but advanced `ref_pos` unclamped, so
it would emit `alignment_end > contig.len()` the moment it became reachable — which Mi3's fix does.
**Applied:** clamps the position, plus a `debug_assert!` on the invariant.

**Mi2 — `summary_array!` read fields by name, so a ninth `RunSummary` field would silently leave the
parity claim.** `RunSummary::merge` destructures exhaustively for exactly this reason. **Applied:**
replaced by a named `SummaryCounters` struct built by two exhaustive destructures — which also
removes the `[u64; 8]` / parallel `SUMMARY_FIELDS` pairing that could name the wrong counter, and the
coverage test's `summary[2]`/`[3]`/`[4]`/`[7]` index access.

**Mi3 — reads never came within 65 bp of a contig end**, so the fold was never differentiated at the
one arithmetic boundary a fetch can get wrong, and four bounds guards were provably inert. Max
`alignment_end` measured 95 on a 160 bp contig. **Applied:** 1-in-8 reads are placed at the far end.

**Mi4 — the `mate_lookup_window` eviction path fired in 0.9% of cases and was not asserted.**
It is the path where a pair silently degrades to two solos. **Applied:** the window is drawn small
enough to bite (1..8 rather than 1..40) and `mate_evictions > 0` is asserted; measured 61 → 221.

**Mi5 — the catch-unwind/collect/summary block appeared four times**, twice inlined in the real-data
test, and the inlined pair had already drifted. **Applied:** one generic `drive` function, used by
all four sites.

**Mi6 — `into_production` had zero runtime coverage** in the default test run; its only caller is the
`#[ignore]`d real-data test. A lossy conversion would surface as "the walkers disagree" — the wrong
diagnosis. **Applied:** `into_production_moves_every_field_to_its_counterpart` and
`the_two_conversions_round_trip_every_field_the_walk_reads`.

**Mi7 — the coverage test hard-coded `CASES_PER_SEED`** while the differential used
`cases_per_seed()`, so a soak widened the comparison but not the coverage floor. **Applied.**

**Mi8 — errors were rendered with `to_string()`**, which is weaker than `{:?}`: an `io::ErrorKind`
inside `WalkerError::Fasta` is invisible through `Display`, and `prepared_read.rs` already uses
`{:?}` for the same distinct-but-identical-types problem. **Applied.**

**Mi9 — `production_reads_len` counted records, not reads**, disagreeing with its name, its parameter
name and its doc, and was a `.len()` pass-through. **Applied:** deleted.

**Mi10 — `ReadGroupId(0)` was written where `PLACEHOLDER_READ_GROUP` exists for exactly this.**
**Applied.**

### Nits

Applied: the module doc's "one `MockFasta`, lent to both" was inaccurate (two are built from the same
bytes; the type is stateless, so it is equivalent — but the real-data test is the one that genuinely
lends one). The "small minority" characterisation of the panic rate is replaced by the measured
figure. Not applied, recorded: `mq_log_err` is drawn independently of `mapq`, so the fold never sees
the two agree; the reference never contains `N`; the sort comment's mechanism is one step longer than
stated (the cap truncates in `ActiveReads` iteration order, which equals admission order only until
the first `swap_remove` expiry — the conclusion holds); there is no `case_at(seed, index)` entry
point, so replaying case *N* means replaying `0..N`.

## 7. Out of scope observations

- **`Cargo.toml` sets `panic = "abort"` on `[profile.release]`**, which `[profile.bench]` inherits —
  the profile `cargo test --release` builds under. Cargo documents the setting as ignored for test
  and bench targets (libtest requires unwinding), so `catch_unwind` still works; and if it did not,
  the process would abort loudly rather than pass silently. Verified from the profile only.
- The reliability pass confirmed **no shared mutable state** in either walker tree
  (`grep` for `thread_local` / `static mut` / `OnceLock` / `LazyLock` returns nothing), so an unwound
  walker leaves nothing behind and later cases are unaffected.

## 8. Missing tests to add now

All added: `both_walkers_report_the_same_error_on_the_same_malformed_input`,
`the_adaptor_filter_changes_the_records_the_walk_emits`,
`into_production_moves_every_field_to_its_counterpart`,
`the_two_conversions_round_trip_every_field_the_walk_reads`, plus the `mate_evictions` floor.

## 9. What's good

- **`catch_unwind` is used correctly**, and the reviewer verified it rather than assuming: `Vec::push`
  is exception-safe so the record prefix survives, `summary` is set last so it stays `None` on a
  panic, and **no panic hook is installed** — which is the right call, since `set_hook` is
  process-global and libtest runs tests on parallel threads.
- **A shorter-but-prefix-equal stream cannot pass** — the length assertion is unconditional and runs
  before the early return.
- **The B2 table is real.** Two rows were re-run independently and reproduced to the seed, the case
  index, the stream item and the record counts.
- **`into_production` is lossless for everything the walk reads**, confirmed by grep: ng's copied
  walker never reads `read_group` on any code path.
- **The generator is a pure function of its seed** — no hashing, no `HashMap` iteration, no clock, no
  thread state; case *N* under seed *S* is invariant to `PVC_PARITY_CASES`.

## 10. Commands to re-verify

`cargo fmt --all -- --check`; `cargo clippy --all-targets --all-features -- -D warnings`;
`cargo test --all-targets --all-features`; `PVC_PARITY_CASES=5000 cargo test --release --lib
ng_walks_identically_to_production`; the four real-data invocations in `parity.rs`'s doc comment.

### Author response convention

Address each finding by identifier with `fixed in <commit>` / `disputed because …` / `deferred to
<issue>` / `won't fix because …`.
