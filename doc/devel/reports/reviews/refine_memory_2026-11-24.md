# Code Review: refine_memory
**Date:** 2026-11-24
**Reviewer:** rust-code-review skill (single reviewer, no category fan-out; see §1)
**Scope:** correctness review of the three refinement-memory commits 790f4a69, 8c9ce22e, d6c543a3
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** three commits on `main` that cut the memory refinement (the
  neighbour-evidence pass) holds per window. The pass holds a window — a core it rewrites plus a
  margin either side read for context — and these commits shrink what each record and each sample
  costs in that window, stream the core's rewritten lines out 1,024 at a time, and, in the last
  commit, read only a sample of the records far from the core ("thinned margins") while holding
  the others compressed.
  - 790f4a69 — f32 sums and posteriors (the largest entry read back as one less the others),
    results made 256 samples at a time, compact phase, compact read alleles, lines counted at
    compressed length, `--refine-window-floor`.
  - 8c9ce22e — likelihoods as a `u16` PL through a lookup table, packed path bits, three-byte read
    alleles, core lines written 1,024 at a time, walk terms recomputed where read.
  - d6c543a3 — token-bucket thinning beyond 50 kb of the core, skipped records zstd-packed, template
    choice over walked records only, `--refine-margin-density`, `--refine-core-floor`, resume
    manifest fields.
- **Reviewed against:** `main` at 822432bb (no later commit touches the reviewed files).
- **In-scope files:** `src/neighbour_evidence/{stage,windows,refine,records,phase,line,model,templates}.rs`,
  `src/read_phase/fragments.rs`, `src/vcf/writer.rs`, `src/cli/{calling_run,resume_calling,call_from_psps}.rs`,
  `src/run/paralog_filter/resume.rs`.
- **Deliberately out of scope:** `examples/refine_replay.rs`, `examples/spill_slice.rs`,
  `examples/vcf_genotype_diff.rs` (measurement tools, not shipped code); the statistical choices
  (density 5 a kb, 50-kb near margin), which the implementation report measured.
- **Categories:** the skill's per-category fan-out in worktrees was **not** run; one reviewer read
  the diffs and the surrounding code with the brief's priorities (correctness, the memory count,
  the input range, resume compatibility) and confirmed the decisive findings with a probe crate
  outside the source tree (§10).

## 2. Verdict

**Request-changes.** One crash (B1) is reachable with the default `--refine-margin-density` on any
run that refines more than one contig; it is a one-line guard. The memory count misses two things
that together plausibly explain most of the 16–22% by which real memory exceeds it (M1, M3), and
thinning added work on the writing thread that grows with every held line and every sample (M2).

## 3. Execution status

- Ran: a probe crate at `tmp/review_2026-11-24_refine_memory/probe/` that depends on the library by
  path (with the repository's `Cargo.lock` and its `[patch.crates-io]` vendored noodles), driving
  `NeighbourEvidenceStage` through its public API. Native `cargo build`, macOS arm64. Outputs quoted
  in the findings.
- Ran: a standalone `rustc` check of how `Option<Vec<T>>` collection sizes its vector
  (`tmp/review_2026-11-24_refine_memory/cap.rs`), quoted in M3.
- Not run: `cargo fmt --check`, `cargo clippy`, `cargo test`, `cargo doc`, `cargo audit` (asked not to
  run the whole suite; the implementation report states 5,363 passed, 2 release-only failures that
  also fail on main).
- Findings labelled "Needs verification": 2 (M1's size estimate, Mi3).

## 4. Open questions and assumptions

1. **Is the evidence-nearby rule meant to see only walked records under thinning?** (Mi1.) The
   spec's thinned-margins paragraph says only template choice and the walks change. I assume the
   rule was meant to read every record, as before.
2. **The margin-thinning report's simulated-cohort table shows "1 a kb, not packed" and "1 a kb,
   packed" writing different genotypes** (outbred wrong 19,506 against 19,553; F2 no call 33,620
   against 35,287). Packing is storage and should change no genotype. I assume the unpacked arm was
   built before template choice was restricted to walked records; if it was not, packing changes the
   output and that needs tracing. Not verified.
3. **Is a dependency bump of zstd allowed to move genotypes?** (Mi3.)

## 5. Top 3 priorities

1. **B1** — guard the empty site span in `compress_held_lines`; the default density panics on the
   second contig when the lines being compressed hold no PASS site.
2. **M1** — count the lines and records the stage keeps to the left of a shrunken window; at the
   margin floor they are about 16% on top of the count in an unthinned window.
3. **M2** — stop rebuilding the joined-record list and the skipped-bytes totals on every pushed line;
   measured 73% more wall time at 10 samples.

## 6. Findings

### Blocker

- `B1: src/neighbour_evidence/stage.rs:613` — **[Blocker]** Compressing a batch with no site panics
  on any contig after the first, with thinning on (the default)
- **Confidence:** High (reproduced).
- **Problem:** `compress_held_lines` asks `site_span(from..held.len())` for the sites among the lines
  it compresses. When none of those lines is a site, `site_span` returns `0..0`, and the next line
  computes `sites.start - first_site`. `first_site` is the run-wide number of the first held site:
  zero on the first contig, and above zero on every later contig and after `forget_what_no_window_reads`
  drops sites. A line is not a site whenever `site_of` refuses it — a non-PASS line, `ALT` = `.`, a
  line without `PL`/`GP`/`AF`. So a short contig whose records are all filtered (an unplaced
  scaffold with a few `LowQual` records), or a contig whose last batch before `finish` is filtered,
  reaches it.
  ```rust
  let sites = self.site_span(from..self.held.len());   // 0..0 when no line is a site
  ...
  self.sites.pack(sites.start - first_site..sites.end - first_site, ...)   // 0 - first_site
  ```
  Probe: 20 PASS sites on contig 0, then one `LowQual` line on contig 1, then `finish`:
  ```
  density 0: 21 lines out
  thread 'main' panicked at src/neighbour_evidence/stage.rs:621:29:
  attempt to subtract with overflow
  ```
  With overflow checks off, as `[profile.release]` has them, the subtraction wraps and the slice
  panics instead (`panic = "abort"` then aborts the run):
  ```
  thread 'main' panicked at src/neighbour_evidence/stage.rs:1127:21:
  range start index 18446744073709551596 out of range for slice of length 0
  ```
- **Why it matters:** a refining run with the default density stops on ordinary input. No test
  reaches it: the two-contig stage test (`run_stage`) runs with `thinning: None`.
- **Suggested fix:** skip packing when the span is empty, and add the two-contig test with thinning
  on and a contig of filtered lines only (§8).
  ```rust
  if self.thinner.is_some() {
      let sites = self.site_span(from..self.held.len());
      if !sites.is_empty() {
          // ... as now
      }
  }
  ```

### Major

- `M1: src/neighbour_evidence/stage.rs:994` — **[Major]** Lines and records held left of a shrunken
  window are not counted
- **Confidence:** High for the mechanism; Medium for the size (Needs verification: measure peak
  against the count with this term added).
- **Problem:** after each window, `forget_what_no_window_reads` keeps every line within
  `sizes.margin_bases` — the *configured* margin, 1 Mb — of the next core, because the next window
  may plan the full margin. When the budget shrinks a window, `plan_window` counts only
  `context_start..context_end`, which starts the *planned* margin before the core. At the floor
  (500-kb margins) the 500 kb between the two are held — their compressed lines, their records'
  likelihoods and shares, their read alleles, and with thinning any records the previous window
  unpacked — and counted by nobody. Lines read past `context_end` while the larger sizes were being
  ruled out are a second, smaller uncounted stretch (at most about two ladder steps, 100 kb at the
  default floor).
- **Why it matters:** the budget is the run's only guard against running out of memory, and this
  gap opens exactly when the budget binds. Using the implementation report's breakdown of the
  failing floor window (per sample per record: 76.5 bytes solve state, 46 read alleles, 24.3
  likelihoods and shares, 8.1 compressed line, 28.4 core results), an unthinned floor window
  counts 1.5 Mb × 154.9 + 0.5 Mb × 28.4 ≈ 246.6 units, and the uncounted 500 kb holds
  0.5 Mb × 78.4 ≈ 39.2 more: **about 16% above the count**, the low end of the measured 16–22%.
  (Assumes records evenly dense across the 2 Mb; thinned windows hold most of that stretch packed,
  so the share there is smaller.)
- **Suggested fix:** count what the stage holds, not only what the window reads. Keep a third
  running total of each line's held cost without solve state (`held` minus
  `samples × held_solve_bytes`), and add its span `0..context_start` and `context_end..held.len()` to
  `plan.bytes`. The extra term only grows as more lines are read, so the planner's "over budget on
  what has been read" rule stays sound. Alternatively, after a window is planned with margin `m`,
  hold nothing older than the largest margin any later window can still choose — but the ladder can
  grow back to 1 Mb, so counting is the simpler correct answer.

- `M2: src/neighbour_evidence/stage.rs:697` — **[Major]** Every pushed line rebuilds two
  whole-window structures on the writing thread when thinning is on
- **Confidence:** High (measured).
- **Problem:** `release` runs once per pushed line, and each iteration of its loop starts with
  `joined_lines()` and `skipped_bytes_besides()`. The first allocates a vector one entry per held
  site and walks every record of every sample's last written block (records within 1 Mb of each
  sample's last phased record); the second clones or rebuilds the held lines' running totals. Before
  this commit a push that did not complete a window cost a few binary searches. These run on the
  thread every line passes through, while the calling threads wait.
  Probe, one contig of 300,000 records 30 bases apart, 10 samples, 4 threads; density 1,000,000 a kb
  walks every record, so the solve's work is the same as unthinned and only this overhead differs:
  ```
  density 0: total 25.2 s
  density 1000000: total 43.7 s
  ```
  The last-block walk grows with the sample count, so at 2,176 samples the overhead per push is
  larger than measured here.
- **Why it matters:** wall time on every refining run with thinning on, growing with records held
  times samples; not yet seen genome-wide (the report says it has not been run genome-wide).
- **Suggested fix:** `last_blocks` changes only when a window is solved, and held lines change only
  by appending or by `forget_what_no_window_reads`. Keep `joined` and the adjusted skipped totals as
  stage fields: rebuild them after `remember_last_blocks` and after `forget_what_no_window_reads`,
  and in `push_writing` append the new line's skipped cost to the adjusted totals as to
  `skipped_bytes` (a new line is never in a last block).

- `M3: src/read_phase/fragments.rs:493` — **[Major]** Three-byte read alleles are held in a vector
  up to twice as long as counted
- **Confidence:** High for the capacity; Medium for its share of the measured excess.
- **Problem:** the common form, `Tight`, is built by collecting an iterator of `Option`s:
  `compact.iter().copied().map(TightReadAllele::of).collect::<Option<Vec<_>>>()`. That adapter
  reports no lower size bound, so `Vec` grows by doubling and keeps the slack. `heap_bytes` counts
  `len × 3`. Checked with the same shape of code:
  ```
  len 1000 cap 1024
  len 1025 cap 2048
  len 5000 cap 8192
  ```
  One vector holds every sample's read alleles at a record (thousands of entries at 2,176 samples),
  so on average the capacity is about 1.44 times the length, and up to 2. In a thinned run the
  records packed and unpacked are rebuilt at their exact length (`read_from` collects from a range),
  so the slack stays on the records the far margins walk and on every record of an unthinned run.
- **Why it matters:** read alleles were 46 of about 183 counted bytes a sample a record in the
  failing window; slack on the entries is several per cent of the window, uncounted, on top of M1.
- **Suggested fix:**
  ```rust
  let tight: Option<Vec<TightReadAllele>> =
      compact.iter().copied().map(TightReadAllele::of).collect();
  let entries = match tight {
      Some(mut tight) => {
          tight.shrink_to_fit();
          PackedEntries::Tight(tight)
      }
      None => PackedEntries::Compact(compact),
  };
  ```
  or check every entry fits first and build `Tight` with `Vec::with_capacity(compact.len())`.

### Minor

- `Mi1: src/neighbour_evidence/refine.rs:1580` — **[Minor]** With thinning, the evidence-nearby rule
  sees only walked records
- **Confidence:** Medium (the behaviour is certain; whether it was intended is open question 1).
- **Problem:** `final_solution` decides whether a sample without reads at a core record is written
  (`has_evidence_nearby`, §4.2: informative records within the sample's segment length `S`, up to
  300 kb, on both sides) from `informative_positions(records, layout, sample)`, which lists only
  `layout.record_of_site` — now only the walked records. Beyond 50 kb of the core, a sample's
  informative records that the bucket skipped no longer count as its evidence. Every record keeps its
  rows (`take_numbers` leaves them), so the full list is available.
- **Why it matters:** a sample with sparse reads — the low-coverage end the project commits to —
  whose nearest informative record on one side lies 50–300 kb past the core edge and was skipped is
  now written `./.` where unthinned it was called. The report's claim that "every record is refined
  as before" does not hold for this rule.
- **Suggested fix:** list positions over every site, as before thinning:
  ```rust
  fn informative_positions(records: &[SiteRecord], sample: usize) -> Vec<u64> {
      records.iter()
          .filter(|record| record.is_site() && record.row(sample) == SampleRow::Informative)
          .map(|record| record.position)
          .collect()
  }
  ```

- `Mi2: src/cli/resume_calling.rs:181` — **[Minor]** Resuming with a new `--refine-window-floor`
  keeps the old core floor
- **Confidence:** High.
- **Problem:** a fresh run's core floor defaults to the window floor
  (`calling_run.rs`: `core_floor.unwrap_or(window_floor)`). On resume, `with_window_floor` changes
  only the margin floor; the core floor read from the manifest (or, from an older manifest, the old
  margin floor) stays. Resuming a default run with `--refine-window-floor 800000` lets cores shrink to
  500 kb under 800-kb margins, which a fresh run with the same option would not, and the header
  gains `core_floor=500000`.
- **Why it matters:** the same options give different windows, and so different genotypes, depending
  on whether the run was resumed.
- **Suggested fix:** when `--refine-window-floor` is given and `--refine-core-floor` is not, set the
  core floor to the new window floor if the manifest's core floor equalled its margin floor (the
  default), or always — and say which in `--help`.

- `Mi3: src/neighbour_evidence/stage.rs:627` — **[Minor]** Window placement, and so genotypes, now
  depend on zstd's compressed sizes
- **Confidence:** Medium (Needs verification: whether the owner accepts this dependency).
- **Problem:** since 790f4a69 a held line counts at its compressed length, and since d6c543a3 a skipped
  record counts at its packed length. Both feed `plan_window`, which chooses window sizes whenever
  the budget binds. The code comments guard against `size_of` moving a genotype; a zstd upgrade
  (`zstd-sys 2.0.16+zstd.1.5.7` in `Cargo.lock`) can change compressed lengths and move a window
  edge the same way. Thread count and platform do not matter: each line is compressed alone, single
  threaded.
- **Why it matters:** the project promises the same bytes from the same input; a dependency bump
  would change output only on budget-bound windows, the case least likely to be in the tests.
- **Suggested fix:** either state the dependency in the spec (§6) and pin zstd's minor version with a
  test that compresses a fixed line and asserts its length, or count lines and packed records at an
  estimate independent of the compressor.

### Nits

- `src/vcf/writer.rs:448` — the doc comment of `write_one_line` was inserted below `place_of`'s, so
  `write_one_line` carries both and `place_of` lost its own. Move the four-line block above
  `place_of`'s doc down to `fn place_of`.

## 7. Out of scope observations

None beyond the findings above.

## 8. Missing tests to add now

- `NeighbourEvidenceStage::compress_held_lines` —
  `a_later_contig_of_filtered_lines_only_passes_through_with_thinning_on`: two contigs, the second
  holding only non-PASS lines, `thinning: MarginThinning::parse("5")`; every line comes out once, in
  order. Catches B1.
- `NeighbourEvidenceStage` — `the_thinned_stage_writes_what_the_unthinned_one_writes_at_a_saturating_density`:
  run `run_stage` at density 1,000,000 and at 0 and assert byte-identical output. Guards the
  packing/unpacking state against drift (every record is walked, so output must match).
- `final_solution` — `a_sample_whose_far_evidence_was_skipped_is_still_written`: a window where one
  sample's only informative record right of the core lies 100 kb out and is not walked; with Mi1
  fixed it is written as unthinned.
- `PackedRecordFragments::from` — `tight_read_alleles_are_held_at_their_length`: assert
  `capacity() == len()` for a record of 1,025 entries. Catches M3.

## 9. What's good

- The planner's decision is made only on fully compressed counts (`release` compresses and re-plans
  before accepting a plan), so window placement does not depend on when a compression batch fell
  ([stage.rs](../../../../src/neighbour_evidence/stage.rs#L712)).
- The thinner is decided once per record as it arrives, in integer units, so the walked set depends
  on neither the window nor the thread count
  ([stage.rs](../../../../src/neighbour_evidence/stage.rs#L1062)).
- Every compact form round-trips at its boundaries and has a test at each: offsets 65,535 and
  65,536, alleles 3/4 and 255/256, Phred 63/64
  ([fragments.rs](../../../../src/read_phase/fragments.rs#L583)); the PL table is checked bit for bit
  over every `u16`.
- New manifest fields are `#[serde(default)]` and read an older manifest as the run it was:
  unthinned, core floor at the margin floor
  ([resume.rs](../../../../src/run/paralog_filter/resume.rs#L361)).

Checked and found correct: the f32 posterior read back with its largest entry as one less the
others (ties to the first, as `score_best_genotype`; all-zero rows are never read); PL clamping at
`u16::MAX`; packed path bits up to 64 paths; the 1,024-line streaming (order kept, each line written
once, a write error stops the run before `core_start` advances, and the writer is then dropped);
join records always walked; `skipped ≤ held` per line, so the planner's monotonicity holds.

## 10. Commands to re-verify

```
cd tmp/review_2026-11-24_refine_memory/probe
CARGO_TARGET_DIR=../target cargo build
../target/debug/refine_probe 0     # 21 lines out
../target/debug/refine_probe 5     # B1: panics at stage.rs:621
CARGO_PROFILE_DEV_OPT_LEVEL=3 CARGO_TARGET_DIR=../target-o3 cargo build --bin timing
RAYON_NUM_THREADS=4 ../target-o3/debug/timing 0 300000 30         # M2
RAYON_NUM_THREADS=4 ../target-o3/debug/timing 1000000 300000 30   # M2
```

### Author response convention

Address each finding by its identifier with `fixed in <commit>` / `disputed because …` /
`deferred to <issue>` / `won't fix because …`. Answer the open questions in §4 first.
