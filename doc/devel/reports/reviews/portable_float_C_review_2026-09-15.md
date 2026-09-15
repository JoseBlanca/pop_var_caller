# Code Review: portable_float_C
**Date:** 2026-09-15
**Reviewer:** rust-code-review skill (orchestrator, single pass — see §1)
**Scope:** uncommitted diff of Milestone C (`doc/devel/implementation_plans/portable_float.md`): the SNP/indel fit's chunks made independent of the pool's width, a cross-platform digest test, a macOS CI job, the identity oracle extended to the fit, two new review checklists on floating-point portability, and the plan's Milestone C/D text
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** `git diff HEAD` in the worktree `pop_var_caller-portable-float` on top of
  `6d30e919` (B7 report), plus the untracked `src/cli/cross_platform_digests.rs`,
  `ai/skills/rust-code-review/code_review/float_portability.md` and
  `ai/skills/rust-performance-review/performance_review/float_portability.md`. 8 tracked files,
  +155 / −35.
- **In-scope files:** `src/parameter_estimation/joint/fit.rs`; `src/cli/cross_platform_digests.rs`
  and `src/cli/mod.rs`; `.github/workflows/ci.yml`; `scripts/promote_ng_oracle.sh`; the two new
  checklists and the edits to both review `SKILL.md` files and `performance_review/hot_loops.md`;
  `doc/devel/implementation_plans/portable_float.md`.
- **Deliberately out of scope:** `scripts/float_exp_table.py` and `tmp/` (a later step), and the
  runs under `tmp/oracle_C/` still in progress.
- **Categories:** the sub-agent fan-out was **not** run, because the dispatcher forbade `cargo` and
  builds while measurements run on this machine, which removes what worktree isolation is for. One
  pass covered `float_portability` (the new category, applied to the diff itself),
  `unsafe_concurrency` (rayon), `refactor_safety`, `tooling` (CI, shell), `extras` (does the diff
  match the plan) and doc accuracy. The performance side of the chunk change was reviewed against
  `rust-performance-review`'s `concurrency` and `float_portability` rules.

## 2. Verdict

**Approve-with-changes.** The fit change is correct. `absorb` treats the three per-position lists
properly on every pass: they are empty unless the pass collects them, so appending an empty list
is a no-op, and `collect()` over an indexed parallel iterator keeps chunk order. No other place in
`parameter_estimation`, `calling`, `paralog` or `run/paralog_filter` joins floating-point totals
in an order the pool chooses (§7). On the four tomato samples at about 3× the fit is no slower,
and it writes one file on macOS and Linux (`2faba12b`, verified here). The digest test's
arguments match the command structs, and nothing it hashes depends on a path, a time or the
machine. A crate version bump also moves the calls checksum, which its documentation omits (Mi5).

Three things to change before committing:

- **The fixed chunk size makes the joint-fit bench's SNP/indel cases run on one or two chunks**
  (M1). Milestone D plans to compare that bench against the Milestone B build, so that comparison
  would measure the lost parallelism rather than the new `exp`.
- **The chunk fix moves the fitted parameters file** (7 of 574 lines on tomato). It is recorded
  nowhere in the plan. The plan's own principle 4 says a change that moves output is its own
  commit (M2).
- **One Blocker rule in the new code-review checklist is wrong** (M3). Taken as written, it files
  a Blocker against `genetics::lgamma`, which is portable.

## 3. Execution status

- **Commands run by the reviewer** (host, read-only):
  - `git diff`, `git status --short`, `git show 75722336`.
  - `md5 -r` on `tmp/chunkfix/fit_{macos,linux}.toml` and `tmp/ab_chunk/{macos,linux}/fit_{float,chunkfix}.toml`:
    both `chunkfix` fits and both `tmp/chunkfix` fits are `2faba12b89f6136283321b0ad65260f0`; both
    `float` (Milestone B) fits are `0fb3e50df25af8a88fc304de108a59b3`. `diff` of the macOS pair: 7
    of 574 lines differ (a read group's error multiplier, three inbreeding coefficients, the two
    genome-wide concentrations, one slippage row).
  - `tmp/ab_chunk/{macos,linux}/ab.tsv`, read in full. macOS, three rounds: B build 416.7, 442.5,
    430.0 s (mean 429.7); chunk fix 423.5, 447.9, 428.8 s (mean 433.4), +0.9% with overlapping
    ranges. Linux, **two rounds only**: B build 653.6, 636.0 s; chunk fix 647.4, 632.8 s, −0.7%.
  - `tmp/ab_chunk_hostload.tsv` (258 samples, one per 20 s), summarised with
    `uv run --no-project python`: load 4.5 to 24.3, median 13.9. Apart from the fits and the Linux
    VM, the busiest process in any sample took at most 1.5 cores of 18 (`mediaanalysisd`); 86
    samples had some other process above half a core.
  - `tmp/c_gates.log`: the container gates and both test suites passed on this tree — Linux 4,734
    library tests + 64 others, macOS 4,733 + 64, 0 failed. Two more than B7's counts on each
    platform, which are the two new tests.
  - Greps over `src/` for `par_iter`, `into_par_iter`, `par_chunks`, `.reduce(`, `.fold(`,
    `current_num_threads`, `ThreadPoolBuilder`, `thread::spawn`, `HashMap`/`AHashMap` (§7).
- **Commands not run:** every `cargo` command and both oracles, by instruction.
- **Needs verification:** the pre-conversion digests `2c36fb95…` (macOS) and `c910c673…` (Linux)
  and the "three different files at one, four and eight threads" measurement have no log under
  `tmp/` (Mi6).

## 4. Open questions and assumptions

1. **How small does a real census get?** The shipped target is 2,000,000 positions
   (`bindings.rs:888`), thinned from the analysed domain (`loci.rs:848`). A run over a small
   targeted panel has fewer. The size of M1's production effect depends on it.
2. **Is the Linux `check` job on GitHub red today?** B6's container run had 3 example-test
   failures and `PROJECT_STATUS.md` records the same three as platform-independent. If CI runs that
   job, `--examples` already fails on `main` (Mi3).

## 5. Top 3 priorities

1. **M2** — commit the chunk fix on its own, with the checksum it moves, and record it in the plan.
2. **M1** — decide Milestone D's bench baseline (the chunk-fix build, not the B build) and whether
   16,384 stays.
3. **M3** — correct the checklist's Blocker rule before any review applies it.

## 6. Findings

### Major

- **M1:** [src/parameter_estimation/joint/fit.rs:1798-1801](../../../../src/parameter_estimation/joint/fit.rs#L1798-L1801) — **[Major]** A census below 16,384 positions now runs on one core, and every SNP/indel case of the joint-fit bench is one or two chunks
- **Confidence:** High on the mechanism; the size is unmeasured.
- **Problem:** Before, a census was cut into `positions ÷ threads` chunks when that was below 16,384,
  so small censuses still used every thread. Now the chunk is always 16,384. The bench
  `benches/ng_joint_fit_perf.rs` runs in a pool of four and draws 5,000 or 20,000 positions
  (lines 504 and 544). Before: four chunks in every case. After: one chunk at 5,000 positions (four
  of the five SNP/indel cases) and two unequal ones (16,384 + 3,616) at 20,000. At 2,000,000
  positions nothing changes: 123 chunks either way up to 122 threads, and the tomato timings above
  agree.
- **Why it matters:**
  - Plan step D2 measures the new `exp` with "the joint-fit … benches" against "the Milestone B
    build". The B build parallelises those cases and this one does not. The SNP/indel cases would
    report a slowdown of up to about the pool width that has nothing to do with `exp`.
  - A production run whose census is under `16,384 × threads` positions loses parallelism in the
    E-step. At a thousand samples, one chunk is a thousand samples' work.
- **Suggested fix:**
  1. For D2, compare against the chunk-fix build, committed first (M2), and re-take that build's
     bench numbers. Say in the D report that B7's bench figures are not comparable for SNP/indel.
  2. Decide the chunk size by measurement, not by the old comment. A smaller fixed size, for
     example 2,048, keeps small censuses parallel. The per-chunk overhead is one binary search per
     sample. The cost is memory: the collected `Statistics` are about 450 bytes per sample per chunk
     (`genotypes` is samples × 16 nodes × 3 `f64`). So 977 chunks at 2,000,000 positions would hold
     440 MB at 1,000 samples. To bound that, keep the fixed chunk size but fold in fixed-size batches
     of chunks, in order:
     ```rust
     const CHUNKS_PER_BATCH: usize = 64;
     let mut into = empty();
     // reserve as now
     for batch in bounds.chunks(CHUNKS_PER_BATCH) {
         let chunks: Vec<Statistics> = batch.par_iter().copied().map(one_chunk).collect();
         for chunk in &chunks {
             into.absorb(chunk);
         }
     }
     ```
     The order of additions is the same at every pool width, and memory is bounded by the batch size.
     Measure the bench at 5,000 and 20,000 positions, and the tomato fit, before choosing.

- **M2:** [doc/devel/implementation_plans/portable_float.md](../../../implementation_plans/portable_float.md) — **[Major]** The chunk fix moves the fit's output and is not in the plan, and it is bundled with Milestone C
- **Confidence:** High
- **Problem:** The fix changes the tomato fit from `0fb3e50d` to `2faba12b`, 7 of 574 lines. B7's
  report and §3.2 record `0fb3e50d` as "the same file on macOS and Linux". Neither C's step list nor
  §8 mentions the fix. §3 principle 4 says: "A step that moves output silently is its own commit."
  Principle 5 says the oracle script is reused "unchanged", and C3 changes it. Neither departure is
  recorded.
- **Why it matters:** the baseline C3 is about to record includes the chunk fix. Someone bisecting
  a moved fit later finds the move inside a commit described as CI, tests and skills. The B7
  figure then disagrees with the oracle without explanation.
- **Suggested fix:**
  - Commit `fit.rs` alone first. The message should give the before and after checksums on both
    platforms, the 7 lines, the interleaved timings (macOS +0.9% over three rounds, Linux −0.7% over
    two), and the pool-width defect it fixes.
  - Add to §8: *"C2 found the SNP/indel fit's output depended on the pool's width. Its chunks were
    sized from the thread count and joined by rayon's `reduce`. Fixed-size chunks, added in order,
    moved the tomato fit from `0fb3e50d` to `2faba12b` (7 of 574 lines), identically on macOS and
    Linux; the B7 fit checksum is superseded. The identity oracle script was extended (C3), so
    principle 5's 'unchanged' no longer holds."*
  - Also in §8: B7's macOS and Linux fits matched while running on 18 and 8 threads. A cross-platform
    match on real data therefore did not show width independence. Only the fixture at 1, 4 and 7
    threads did.

- **M3:** [ai/skills/rust-code-review/code_review/float_portability.md:14](../../../../ai/skills/rust-code-review/code_review/float_portability.md) — **[Major]** The Blocker rule counts a direct `libm` call as reaching the platform library; `genetics::lgamma` does exactly that and is portable
- **Confidence:** High
- **Problem:** The rule lists "`libm` functions called directly instead of through `crate::float`"
  among "a call that reaches the platform library another way", at Blocker. `libm` is pure Rust and
  is the same on every target: the B6 review said so, and `src/float.rs:22` says `lgamma` "has always
  called libm". `src/genetics.rs:57` calls `libm::lgamma(x)`. A review applying the checklist would
  file a Blocker against correct shipped code. Because the rule is labelled *required*, a reviewer
  may not question it.
- **Suggested fix:** take `libm` out of the Blocker list. If routing matters, add a separate,
  lower rule:
  ```markdown
  - **Nit:** a `libm` function called directly outside `src/float.rs` and `src/genetics.rs`.
    It is portable, but it bypasses the one module where the pinned-bits tests live; route it
    through `crate::float`.
  ```

### Minor

- **Mi1:** [src/parameter_estimation/joint/fit.rs:1759-1762](../../../../src/parameter_estimation/joint/fit.rs#L1759-L1762) — **[Minor]** `expectation_pass`'s doc still describes the streaming sum it no longer has
- **Confidence:** High
- **Problem:** "the chunks have to be held until they can be joined in position order rather than
  summed as they arrive — so the iterating passes use the streaming sum and pay neither." Every
  pass now holds its chunks.
- **Suggested fix:**
  ```rust
  /// **Kept only on the last pass of a run.** Keeping it costs one four-byte value a position.
  /// Every pass, this one included, holds its chunks until all have finished and adds them in
  /// position order, so the totals do not depend on the pool's width.
  ```

- **Mi2:** [src/cli/cross_platform_digests.rs:204-216](../../../../src/cli/cross_platform_digests.rs#L204-L216) — **[Minor]** The pool-width test cannot see a width-dependent join any more, because its census is one chunk
- **Confidence:** High
- **Problem:** The fixture's contig is 600 bases, so its census is under 16,384 positions and the
  E-step is one chunk at every width. The test would catch someone putting back chunk sizes
  derived from the thread count (600 ÷ 4 = 150). It would not catch `reduce` returning over fixed
  chunks, the other half of the defect, because one chunk has nothing to join.
- **Suggested fix:** add a unit test in `fit.rs` whose census spans at least three chunks. The
  bench's cohort drawer already makes 20,000-position cohorts. Fit one of about 50,000 positions
  inside pools of 1 and 4 threads and assert the fitted parameters are equal with `to_bits`. If
  that is too slow for the suite, make the chunk size a field of `JointFitConfig` with 16,384 as
  its default, and use a small value in the test.

- **Mi3:** [.github/workflows/ci.yml:90-93](../../../../.github/workflows/ci.yml) — **[Minor]** The comment blames macOS for three example failures that fail on Linux too
- **Confidence:** Medium (the Linux evidence is the arm64 container, not the x86_64 runner)
- **Problem:** The comment says the three tests "fail on macOS on `main` as well". The B6 review's
  container run had "3 failed that fail identically before the plan". `PROJECT_STATUS.md:1545-1558`
  describes the same three as behaviour and fixture defects in the examples, with no platform
  involved. So the reason to leave `--examples` out is not macOS. By the same evidence, the Linux
  job's `--examples` step fails on `main`.
- **Suggested fix:** reword the comment: *"three example tests (…) fail on every platform, on `main`
  as well (`PROJECT_STATUS.md`, locus dumps). They are research tools outside the portability
  guarantee."* Better, mark the three `#[ignore = "fails on main; see PROJECT_STATUS.md"]` and run
  `--examples` in both jobs.

- **Mi4:** [scripts/promote_ng_oracle.sh:82-88, 148-161](../../../../scripts/promote_ng_oracle.sh) — **[Minor]** A missing baseline reads as seven mismatches, and a failed fit exits with no message
- **Confidence:** High
- **Problem:**
  - `scripts/promote_ng_oracle.baseline` does not exist yet, and the default compares against it.
    `grep -v '^#' "$baseline"` then prints "No such file" and exits 2. The pipeline's status is the
    last `grep -qx`, which is 1. So every one of the 7 lines prints "DIFFERS from …", as if the
    output had moved.
  - Under `set -e`, a failing `estimate-parameters` or `call-from-psps` ends the script right after
    "=== the fit ===". The only explanation is in `fit.log` or `from_fit.log`.
- **Suggested fix:**
  ```sh
  "$bin" estimate-parameters … > "$out/fit.log" 2>&1 \
    || { echo "the fit failed; see $out/fit.log" >&2; exit 1; }
  …
  if [ "$baseline" != none ]; then
    if [ ! -r "$baseline" ]; then
      echo "no baseline at $baseline; set PROMOTE_NG_BASELINE=none to only record" >&2
      exit 1
    fi
    …
      if ! grep -v '^#' "$baseline" | tr ' ' '|' | grep -qxF "$line"; then
  ```

- **Mi5:** [src/cli/cross_platform_digests.rs:28-30](../../../../src/cli/cross_platform_digests.rs#L28-L30) — **[Minor]** A new crate version moves the calls checksum too, and the doc says it moves only the fit's
- **Confidence:** High
- **Problem:** The VCF header carries `##source=ng <version>` (`src/vcf/header.rs:143`), and
  `comparable_vcf` strips only `##commandline` and `##reference`. So a version bump moves
  `CALLS_MD5` as well as the census digest in the parameters file. Separately, the module doc
  (line 36) and the `fit.rs` comment say the defect showed "at one, four and eight threads", while
  the test runs one, four and seven.
- **Suggested fix:** "**A new crate version**, which is written into the parameters file's census
  digests and into the VCF's `##source` line. Re-record both after checking that only those lines
  moved." Say "one, four and eight threads when measured; the test uses seven", or measure at seven.

- **Mi6:** [src/cli/cross_platform_digests.rs:61-66](../../../../src/cli/cross_platform_digests.rs#L61-L66) — **[Minor]** The measurements that give the test its teeth have no log — **Needs verification**
- **Confidence:** Medium
- **Problem:** The pre-conversion checksums `2c36fb95…` and `c910c673…`, the unchanged `CALLS_MD5`
  on that tree, and the three files at 1, 4 and 8 threads are cited in code. None appears under
  `tmp/` (searched for the two checksums and for the module name: only `tmp/ab2/chunk_fix.patch`
  matches). The plan's verification table makes "failing across platforms on the tree before the
  conversion" a Milestone C criterion.
- **Suggested fix:** keep the two test outputs and the width run in, for example,
  `tmp/digests_C/{macos,linux}_c6a4394b.log` and cite that path in the C report.

- **Mi7:** [ai/skills/rust-code-review/code_review/float_portability.md](../../../../ai/skills/rust-code-review/code_review/float_portability.md), [performance_review/float_portability.md](../../../../ai/skills/rust-performance-review/performance_review/float_portability.md) — **[Minor]** Five factual claims are wrong or wider than what was measured
- **Confidence:** High
- **Problem**, each checked against the source:
  1. Code checklist line 5: "one unit in the last binary place of a `ln` flipped a tie". Commit
     `75722336` says **two** units (`…dd29` against `…dd27`), in a table computed with `powf` and
     `ln`.
  2. Code line 13: "glibc, Apple's library and glibc's x86_64 variants … disagree on about one
     argument in 540 to 2,700 for `exp` and `powf`". A2 measured macOS against arm64 Linux. No
     x86_64 variant was measured (B7 §5).
  3. Code line 27: "over 160 regions, macOS and Linux wrote identical VCFs even with the platform
     library in place". That was with **default parameters** (B7 §3.1). Each platform calling with its
     own fit gave 46 differing records of 6,735, one in QUAL (A3, quoted in B7 §3.3). Without the
     qualifier the sentence undersells the risk.
  4. Performance line 14: "libm's `exp` takes 1.3 to 2.1 times as long per call as glibc's". A2's
     Linux rows are 1.26, 2.12 and 0.99, the last past underflow. Say "except past underflow, where
     they are level".
  5. Performance line 18: "the joint-fit bench's draws turned out identical". B7 checked that on
     macOS only.
- **Suggested fix:** correct the five sentences as above.

- **Mi8:** [performance_review/float_portability.md:13](../../../../ai/skills/rust-performance-review/performance_review/float_portability.md), [code_review/float_portability.md:24](../../../../ai/skills/rust-code-review/code_review/float_portability.md), [hot_loops.md:21-22](../../../../ai/skills/rust-performance-review/performance_review/hot_loops.md) — **[Minor]** "a fixed-size chunked (pairwise) sum" is offered as the portable way to vectorise, but a chunked sum with a serial inner loop does not vectorise
- **Confidence:** Medium-high
- **Problem:** Three different things are treated as one:
  - A chunk summed left to right is still one serial chain, so the compiler cannot vectorise it.
    The fit's 16,384-position chunks buy thread parallelism, not SIMD.
  - What vectorises, and stays portable, is a **fixed number of lane accumulators** written in the
    source: `[f64; 8]` updated over `chunks_exact(8)` and combined in a fixed order at the end.
  - Pairwise summation is about accuracy.

  In `hot_loops.md`, the new "Superseded" bullet sits directly above the unedited text that calls
  `algebraic_*` "the *first* fix to reach for here". A reader following the order meets both.
- **Suggested fix:** replace the checklists' sentence with *"The portable way to get a vectorisable
  float sum is a fixed number of lane accumulators written in the source (e.g. `[f64; 8]` over
  `chunks_exact(8)`, combined in a fixed order), whose order does not depend on the target."* In
  `hot_loops.md`, edit the `algebraic_*` bullet itself: "only where the result never reaches
  output". Do not leave the contradiction marked as kept.

### Nits

- `src/cli/mod.rs:14`: `mod cross_platform_digests;` is private where the model module is
  `pub mod mode_equivalence;`. Both compile to an empty module outside tests. Private is arguably
  better, because rustdoc does not publish a doc for a module with nothing in it. Pick one
  convention.
- `cross_platform_digests.rs:183-186`: `called.contains("PARALOG_POST=")` searches the header as
  well. The header happens to spell it `PARALOG_POST,`. Search the record lines only.
- `promote_ng_oracle.sh:67`: `PROMOTE_NG_FIT=no` or `false` runs the fit. Say "`0` skips it" in the
  doc, or accept those values too.
- `promote_ng_oracle.sh:69-74`: the binary search copies the wrapped oracle's, including picking
  a newer macOS binary inside the Linux container, where `-x` passes and exec fails.
  `tmp/c_gates.sh` works around this with `touch`. It is pre-existing and consistent. A check that
  the candidate runs (`"$candidate" --version >/dev/null 2>&1`) in both scripts would remove the
  trap.
- `promote_ng_oracle.sh:152`: the baseline comparison skips the `#` lines, so a run over a
  different `PROMOTE_NG_REGIONS` reports "DIFFERS" for every checksum, not "different slice". The
  script's header promises runs "cannot be compared by accident". Compare `# regions kept` and
  `# alignments` too.
- `promote_ng_oracle.sh:34-35`: "about 6 on the macOS host and 8 in the Linux container". The
  measured fits are 6.3 min (B7, macOS), 7.0–7.5 min (this change, loaded host) and 9.7–10.9 min
  on Linux. Say "about 7 and 10".
- `ci.yml`: `macos-latest` is arm64 on GitHub's current images, and Intel images are being
  retired. Pinning `macos-15` would make the job name's "(arm64)" a fact rather than a current
  default. Consider `timeout-minutes`.
- Code checklist line 20: `fold_chunks` has a fixed chunk size and yields per-chunk results in
  order. It is deterministic when collected; the defect is a `.sum()`/`.reduce()` after it.
- Code checklist line 23: std's `RandomState` differs per map instance as well as per process.
- Performance checklist: `concurrency.md` has a "Note on overlap with the correctness review";
  this category would benefit from the same, since its rules produce rejections rather than
  Hot-path/Likely findings. Two facts worth a line each:
  - `target-cpu`, LTO and `codegen-units` do not change scalar float results in Rust, because rustc
    never contracts `*`/`+` into FMA. Only `algebraic_*` code and target-dependent chunking are
    exposed.
  - On x86_64 without the FMA feature, `mul_add` becomes a call to the C library's `fma`, which is
    exact but slow.

## 7. Out of scope observations

**Every floating-point parallel site outside `fit.rs`, and whether its result depends on the pool.**
None do:

| site | what runs in parallel | how it is joined | reaches output |
|---|---|---|---|
| `parameter_estimation/joint/ssr_fit.rs:1144` | one climb per (stratum, starting point) | `collect`, then `fold(the_better_walk)` in starting-point order | yes, width-independent |
| `ssr_fit.rs:2167-2177` | `ln_tract` per tract | `collect` into `Vec<f64>`, then serial `sum` | yes, width-independent |
| `ssr_fit.rs:2206`, `2271` | tract likelihood tables; Dirichlet points by row | `collect`; `par_chunks_mut` writes disjoint rows | yes, no summation |
| `parameter_estimation/joint/contamination.rs:621` | one search per library | `collect`; the inner `fold` is over integers | yes, width-independent |
| `cli/estimate_contamination.rs:440` | one census walk per sample | `collect` in input order | yes, width-independent |
| `run/paralog_filter/pass_two.rs:304` | scoring a batch of spilled records | `collect_into_vec`, pushed in order | yes, width-independent |
| `run/cohort_merge/{rounds,parallel}.rs`, `observation_cache.rs` | regions a round, samples | index-ordered; region width follows the pool, and `call_from_psps.rs:171` and the existing serial-cover tests pin calls identical at every width | yes, pinned by existing tests |
| `reference_info.rs:662` | contig MD5 split | integer digest | no float |

- **`HashMap` feeding float totals:** none. The only hash map on a calling path,
  `run/cohort_merge/build.rs:518`, is a lookup that nothing iterates.
- **`fold` over floats:** every `fold` in these modules is a `max`/`min` (order-independent for
  non-NaN values), an integer count or a digest.

**Memory of the collected chunks.** Each `Statistics` is about 450 bytes per sample (16 quadrature
nodes × 3 `f64` in `genotypes`, plus five per-sample totals) and 144 bytes per read group. Holding
123 of them at 2,000,000 positions costs 3.5 MB at 63 samples, 55 MB at 1,000 and 165 MB at 3,000.
They are now held on every pass. The final pass already held them, together with lists of 12 bytes
a position a sample, so the run's peak does not rise. No action unless M1's fix raises the chunk
count.

## 8. Missing tests to add now

- The multi-chunk pool-width unit test in Mi2.

## 9. What's good

- **The fit fix is the smallest correct one.** It deletes the thread-derived size and the `reduce`
  arm, and it keeps `absorb`, the reservations and the output order unchanged. `collect_genotype_posterior`
  implies `collect_noisy_posterior` (line 1806), so reserving the genotype lists outside the noisy
  branch changes nothing.
- **The digest test measures the right artefact.** Its documentation states why a calls-only
  checksum would pass with the defect in place, and cites B7 §3.1 for it. The second test runs each
  width on a fresh fixture in a different temporary directory. When it passes, that also proves no
  temporary path reaches either artefact beyond the two stripped header lines.
- **The test is isolated from the rest of the suite.** `threads: 0` makes `run_call_from_psps` skip
  `build_global` (`call_from_psps.rs:402`). Every parallel site in the walk, fit and call asks
  `rayon::current_num_threads()` inside `pool.install`, so tests that set the global pool cannot
  change this one's width.
- **Sample order is fixed on every path the test takes.** Calling names `one` then `two`
  explicitly. `estimate-parameters` expands the directory with `psps_named`, which sorts by path
  (`psp_inputs.rs:78`).
- **Nothing machine-dependent is hashed.** The parameters file holds `reference_digest` and census
  digests of contents, not paths. `##parametersFile` is a bare file name. `##paralogFilter` holds
  numbers only. Rust's float formatting does not depend on locale.
- **The oracle script's comparison is correct under POSIX `sh` and `set -eu`.** The two-space
  checksum lines become `hash||name`, with no spaces or glob characters left to split. `grep`
  statuses inside `if !` do not trigger `set -e`. The fit's psps follow the `#CHROM` order, as the
  wrapped oracle's own `call-from-psps` does. With `PROMOTE_NG_FIT=0` only the five lines written
  are compared.
- **CI:** the YAML is well formed. The workflow-level `RUSTFLAGS: -D warnings` applies to the new
  job. `rust-toolchain.toml` pins 1.98.0 for it as for the Linux job. `macos-latest` is arm64.
- **Plan text:** the Milestone C amendments and Milestone D match B7. Checked: default-parameter
  calls identical before the conversion (§3.1); `exp` about 20% of the unchanged fit's CPU (A3
  line 147); libm's `exp` 1.3–2.1× glibc's outside underflow (A2). The skills' performance figures
  (fit +27%/+29%, `call-from-psps` +1.7–2.0%, joint-fit bench +7 to +46%, 496 s against 280–288 s,
  `generate-psps` +0.8%/−1.0%) all match B7.
- **The skill wiring follows the existing structure.** Both new files have Purpose, Triggers, Skip
  when and Rules, like `extras.md`. Each `SKILL.md` gains one principle and one triage row marked
  required from 2026-09-15, which is what the owner asked for.

## 10. Commands to re-verify

- The gates, unchanged: `cargo fmt --all -- --check`,
  `cargo clippy --all-targets --all-features -- -D warnings`,
  `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features`,
  `cargo test --lib --tests --all-features` on both platforms.
- `cargo test --lib cross_platform_digests` on both platforms; on the first CI run, the Linux job is
  the first x86_64 evidence. If it fails there, the module doc's rule applies: do not re-record.
- After Mi2: the new unit test at 1 and 4 threads, then the same test with `.reduce` put back over
  fixed chunks, which must fail.
- For M1: `cargo bench --bench ng_joint_fit_perf --features bench-fixtures` on the B build and the
  chunk-fix build, SNP/indel group only, four threads.
- `sh -n scripts/promote_ng_oracle.sh`; then a run with `PROMOTE_NG_BASELINE` pointing at a
  missing file, which must say so (Mi4).

## Fixes applied

| finding | resolution |
|---|---|
| M1 small censuses on one core | measured and redesigned. Fixed 2,048-position chunks in batches of 64 made the bench's SNP/indel cases 47–56% slower; fixed 256-position chunks in a join tree made one tomato fit 28% slower and 99 MiB larger. Chosen: the chunk size from the census length only (128 chunks, at least 256 positions), iterating passes joined in a fixed `rayon::join` tree, the reported pass collected in order. Bench SNP/indel cases +0.5 to +1.7% against `cd056a3f` (macOS, same session); tomato fit +2.2 to +2.8% (macOS, three rounds) and −1.8, −5.7% (Linux), interleaved. D's baseline will be this commit |
| M2 output-moving fix bundled | committed alone as `f04f7c43` with before and after checksums (`0fb3e50d` → `9e0e7a4a`, both platforms, 7 of 574 lines) and timings; plan §8 records it, the superseded B7 checksum, and why B7's 18- and 8-thread match did not show width independence |
| M3 direct `libm` as Blocker | removed; a direct `libm` call outside `float.rs`/`genetics.rs` is now a Nit for routing |
| Mi1 stale doc | reworded |
| Mi2 one-chunk fixture | `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` added in `fit.rs`; with the old thread-sized chunks and `reduce` put back it fails (checked at the 2,048 and 256 chunk sizes) |
| Mi3 CI comment | reworded: the three failures are on every platform measured, and the Linux `--examples` step is likely failing on `main`; runner pinned to `macos-15`, `timeout-minutes: 60` |
| Mi4 oracle robustness | unreadable baseline exits with a message; the fit and the call with the fit exit with a message naming their log; `grep -qxF`; the chosen binary must run; the slice walked is compared first; `PROMOTE_NG_*` must be set inside the container command, which the script's usage now says (found when a run inside the container ignored `PROMOTE_NG_BASELINE`) |
| Mi5 version bump, thread counts | doc says a version bump moves both checksums; the 1/4/8 measurement is cited with its log and the test's seven stated |
| Mi6 unlogged evidence | `tmp/digests_C/{macos,linux}_c6a4394b_final.log` (pre-conversion checksums with the final fit code applied) and `tmp/digests_C/pool_width_cd056a3f.log` |
| Mi7 five claims | all five corrected in the two checklists |
| Mi8 vectorisation | both checklists name fixed lane accumulators as the portable vectorisable form; `hot_loops.md`'s `algebraic_*` bullet itself now restricts it to values that never reach output, and the chunked-sum sentence says a chunk summed left to right does not vectorise |
| nits | `PARALOG_POST=` searched in records only; fit time "about 7 and 10" minutes; `PROMOTE_NG_FIT` doc says exactly `0`; `fold_chunks` and `RandomState` wording; overlap note and the `target-cpu`/`mul_add` facts added to the performance checklist. The module stays private (`mod cross_platform_digests;`), since it holds only tests |
