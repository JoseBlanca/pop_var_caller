# Code Review: neighbour_evidence_phasing
**Date:** 2026-11-02
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** branch `neighbour-pass`, commits 20d7df26..1fbe192d — phased templates for the neighbour-evidence pass
**Status:** Approve-with-changes (the changes are applied: 0919d095 and the commit after it)

---

## 1. Scope

- **Reviewed:** the diff 20d7df26..1fbe192d against its spec,
  [neighbour_evidence_phasing.md](../../ng/spec/neighbour_evidence_phasing.md).
- **In-scope files:** `src/neighbour_evidence/{random,model,templates,refine,stage,fit,mod}.rs`, their
  tests, `src/neighbour_evidence/testdata/phased_window.txt` and its generator
  `doc/devel/ng/research/neighbour_evidence/phaser_reference.py`.
- **Out of scope:** everything outside `src/neighbour_evidence` (pre-existing), the other research
  scripts.
- **Categories dispatched** (one sub-agent each, each in its own git worktree, detached at
  1fbe192d): reliability; errors; float_portability; unsafe_concurrency (rayon in the window solve);
  refactor_safety; module_structure; naming with defaults; idiomatic with smells; extras (diff
  against the spec, the diff's own numbers, hot path). `tooling` was skipped: no `Cargo.toml` change.
  The per-category files are in `tmp/review_2026-11-02_phased-templates/`.

## 2. Verdict

Approve-with-changes. No defect reachable from the pass's real inputs was found; the concurrency
and float-portability reviews found the output independent of thread count and of the maths
library. The findings were test gaps (one Blocker: nothing above two alleles), memory and cost far
above what the spec stated, one reproducibility defect (draws keyed by the run-wide site number),
and stale prose and numbers.

## 3. Execution status

- `cargo test --release` (lib): 5,046 passed, 2 failed — `#[should_panic]` tests on debug
  assertions, failing because the run was release; integration tests pass.
- `cargo clippy --release --all-targets`: clean apart from vendored `noodles-cram`. With
  `--all-features -- -D warnings`: a pre-existing error in test `ng_calling_loop_allocation`.
- `cargo doc --no-deps`: two broken links (fixed, M9).
- `cargo fmt --check`: diffs only in `src/locus_generation/` (out of scope).
- `cargo audit`: not installed.
- Mutations run by the reviewers: reliability 31 (6 survived, 5 of them findings), refactor safety 5
  (1 survived), concurrency 2 (both killed), errors 3 (all changed no behaviour), float portability 1
  (changed bits, no test failed).

## 4. Open questions and assumptions

None left open: each finding below is fixed or disputed with its reason.

## 5. Top 3 priorities

1. **B1** — nothing tested phasing beyond two alleles. Fixed.
2. **M1/M2** — memory per round about 7 times the spec's; the haplotype choice at about 30 ns a
   comparison against the 0.3 the spec assumed. Fixed (double-buffered store written in place; table
   lookups); the cost is not re-measured at scale.
3. **M3** — draws keyed by the run-wide site number, so a contig's draws changed with what was called
   before it. Fixed (contig and position).

## 6. Findings

Each with its disposition. Fixes are in 0919d095 (the code) and the commit after it (tests and
fixture).

### Blocker

- **B1: `templates.rs` / `model.rs` — phasing never tested above two alleles.** Categories:
  reliability, refactor_safety. **Fixed:** the draw, the re-estimation and the haplotype choice at
  three alleles, against hand values (`a_drawn_path_follows_the_only_template_that_fits_at_three_alleles`,
  `the_re_estimation_at_three_alleles`, `the_haplotype_choice_reads_every_allele`). The Python
  reference stays biallelic; the spec says so (§8.1).

### Major

- **M1: `refine.rs` — a round held a second copy of every sample's haplotypes and core quantities**
  (80–110 bytes a sample-record) before merging. Categories: unsafe_concurrency, extras, idiomatic.
  **Fixed:** `HaplotypeStore`, double-buffered; each sample writes its own row from its own thread;
  sums and previous posteriors updated inside the per-sample closure.
- **M2: `refine.rs` — two binary searches per comparison in the haplotype choice**, about 30 ns each,
  500 core-hours at 3,000 samples by extrapolation. Categories: extras, idiomatic, unsafe_concurrency.
  **Fixed:** `WindowLayout`'s record-to-site table; previous posteriors passed aligned with the
  compared records. Not re-measured.
- **M3: `stage.rs` — draws keyed by the run-wide site number.** Category: float_portability
  (cross-category). **Fixed:** `window_key(contig, first core position)`. The spec is corrected: the
  genotypes can still depend on what was called before, through the parameter fit's first windows,
  which predates phasing.
- **M4: `refine.rs` — three index spaces (record, site, core site) as bare `usize`;** a mix-up
  survived the suite because every test window had only sites. Category: refactor_safety.
  **Fixed:** `WindowLayout` holds every conversion; the reviewer's test
  (`a_non_site_record_changes_no_other_phased_record`) is committed.
- **M5: re-fitting in every round, and the fill check fed the posterior, both went unnoticed.**
  Category: reliability. **Fixed:** `only_round_one_fits`,
  `the_fill_check_compares_the_neighbours_prediction`; each was checked to fail under its mutation.
- **M6: the 50-record cap and the stage's window key untested.** Category: reliability. **Fixed:**
  `hidden_records_respect_both_limits` (50 of 3,000; 20 of 400), `window_keys_differ_by_contig_and_position`.
- **M7: the Python fixture too small** (8 samples: the 16-slot cut and the one-a-sample rule never
  compared; hidden lines ignored). Category: reliability, extras. **Fixed:** 20 samples, 500 records;
  hidden records asserted; Python and numpy versions in the fixture header.
- **M8: four `expect`s without the `// PANIC-FREE:` note.** Category: errors. **Fixed.**
- **M9: nothing at run time stated the hidden-record limits; the fill check's docs and log described
  the unphased check.** Categories: naming_defaults, extras. **Fixed:** the header gains
  `hidden_per_window=50,hidden_share=0.05`; `FillCheck`, the module docs and the log say the check
  runs at hidden records with phasing; both `cargo doc` links fixed.

### Minor

- `(Vec<f64>, Vec<f64>, f64)` tuples for posterior, prediction, background — **fixed** (`CoreSite`).
- `Phasing` could hold zero rounds and its "is phased" test was repeated in five places — **fixed**
  (private fields, `new` panics on zero, `is_on`/`averages`/`averaged_rounds`; `phasing_refuses_zero_rounds`).
- New items `pub` though used only inside the pass — **fixed** for `random`, `walk_and_draw`,
  `HaplotypePair`, `compared_records`, `choose_haplotypes`.
- The fit counted hidden records as informative — **fixed**.
- `walk_and_draw` did not check ascending positions; its `# Panics` omitted the zero-weight case —
  **fixed**.
- Phased tests below 30 samples — **fixed** (`a_phased_window_of_a_few_samples_runs`, 2 and 3).
- No cross-platform pin of the pass's bits — **fixed** for one window
  (`a_phased_windows_posteriors_are_the_same_bits_everywhere`, recorded on macOS arm64; Linux must
  confirm it).
- Stale numbers: "867 against 2,795" is 1,015 against 2,795 since up to 50 records are hidden —
  **fixed**; db642c41's "differs by 0.42" is 0.55–0.63 on today's fixture — commit message, left.
- Spec: round 1 is today's pass except at hidden records; hidden reads reach other samples' round-1
  haplotypes through `GP` (masked: no systematic effect); measured cost; the built store's layout —
  **fixed** in the spec.
- `HIDDEN_RECORDS_PER_SAMPLE` renamed `..._PER_WINDOW`; `hidden_records` → `draw_hidden_records`;
  `previous` → `previous_posteriors` — **fixed**.
- **Disputed:** one constant for round-1 template samples and haplotype templates — both are 16 by
  the owner's ruling; one parameter is simpler, and the header prints both.
- **Disputed:** splitting `refine.rs` into `refine/phasing.rs` — the round loop is round 1's solve;
  `WindowLayout` and `HaplotypeStore` now carry the structure; revisit if the file grows.
- **Disputed:** one `transition()` shared by the forward walk, the backward walk and the draw — the
  three agree term by term (idiomatic reviewer), the Python reference catches every transition
  mutation tried (four), and refactoring the forward walk risks the 1e-9 parity for no behaviour.
- **Deferred:** the `A − 1`, `f32` store (spec §5); measuring the cost at 200–3,000 samples.

### Nits

Not enumerated; the ones fixed: `frozen_at` without struct-update syntax; `FillCheck` field docs;
the module index in `mod.rs`; the `SEED` note tying it to the Python reference.

## 7. Out of scope observations

- `cargo clippy --all-features`: test `ng_calling_loop_allocation` does not compile (missing field
  `low_quality_reads`) — pre-existing.
- The SplitMix64 finaliser exists in three other places in the crate; any shared version must stay
  bit-identical to `phaser_reference.py`.

## 8. Missing tests to add now

All the reviewers' proposals with merit are committed (§6). Not added: the 1/4/30-thread test (the
concurrency reviewer checked 1, 3 and 8 on a multi-window stage); `rounds == 2` writing round 2
alone (covered by `only_round_one_fits` and the averaging tests).

## 9. What's good

- The round barrier is enforced by the borrow checker, not by convention: the store being read is
  borrowed by the parallel map and replaced only after it (`refine.rs`).
- The exact comparison with an independent Python reference through a shared random stream turns
  "statistically similar" into "equal to 1e-9", and caught 16 of the reliability reviewer's 25 killed
  mutations alone.
- One round is the old pass by construction, so every pre-phasing test checks round 1.

## 10. Commands to re-verify

- `cargo test --release --lib neighbour_evidence` (60 tests).
- `cargo test --release --lib a_phased_windows_posteriors_are_the_same_bits_everywhere` on Linux
  x86-64.
- `uv run --with numpy python doc/devel/ng/research/neighbour_evidence/phaser_reference.py` regenerates
  the fixture (bit-identical only on the recorded Python/numpy versions).
