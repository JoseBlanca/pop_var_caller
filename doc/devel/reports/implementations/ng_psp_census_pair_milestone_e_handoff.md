# Handoff — Checkpoint D is done; Milestone E is the last one

**Written 2026-09-10**, on branch `census-vs-psp-perf`, in the worktree
`/Users/jose/devel/pop_var_caller-census-vs-psp`. **The tree is clean at `0f4b0646`.** Milestones A
to D are complete, each step implemented, reviewed, mutated and committed. This exists so a fresh
session can pick the plan up without re-deriving anything.

---

## Where the plan stands

The plan is [`psp_census_pair.md`](../../ng/impl_plan/psp_census_pair.md); the design is spec
[`psp_census_pair.md`](../../ng/spec/psp_census_pair.md). Milestone D's commits:

| commit | what |
|---|---|
| `34359259` | D1 — `generate-census` becomes `regenerate-census`, which replaces each psp's trailer |
| `03e5e52f` | D1's review: the blocker (no test could fail if the write never happened), two more inert tests, nine falsified doc claims |
| `c55c55d5` | D2 — a psp whose census this run would write again is skipped, and the cohort is opened before the reference is read |
| `c77fc87b` | D3+D4 — a stopped run does only what is left; a rebuilt psp is the walked one byte for byte |
| `0f4b0646` | `PROJECT_STATUS` at Checkpoint D |

**The reports are** [`ng_psp_census_pair_milestone_d_2026-09-10.md`](ng_psp_census_pair_milestone_d_2026-09-10.md)
and [`../reviews/psp_census_pair_milestone_d_2026-09-10.md`](../reviews/psp_census_pair_milestone_d_2026-09-10.md);
Milestone E appends its sections to both.

---

## What Milestone E's four steps actually need, measured on this tree

Read the plan's own E text first. What follows is what this tree says, which differs from it in
three places.

**E1 — scripts.** Two name the command that no longer exists:

- `scripts/ng_fit_stage_end_to_end.sh` — **broken, not merely stale**: it invokes `generate-census`
  with `--output-dir`, passes `--census` to the fit, and diffs census files. It becomes the
  copy-regenerate-`cmp` of D4, and the shape to copy is
  `a_regenerated_psp_is_the_walked_one_byte_for_byte` in `regenerate_census/tests.rs` — **including
  the part that matters: the copy's trailer has to be emptied first**, or the command skips it and
  `cmp` passes without a rebuild (spec §8; the plan's D4 line records this).
- `scripts/ng_census_route_cost.sh` — names `generate-census`; check what it measures still exists.

**The plan's E1 also names `examples/ng_census_route_cost.rs`'s `write_psp(path, Some(census))`.
That is already done** — the example calls `gatherer.write_psp(&psp_path)` with one argument.

E1's own oracle is the end-to-end script over the six tomato accessions on the two 100 kb intervals:
the same parameters file and the same VCF as the last recorded run.

**E2 — the sidecar's machinery deleted.** `PileupIdentity`, `freshness`, `freshness_by_header`,
`psp_beside`, `census_path_for`, `CENSUS_FILE_EXTENSION`, and what is left of `census_cohort.rs`.
Two things D1's review measured and recorded:

- `CENSUS_FILE_EXTENSION` and `census_path_for` (both in `generate_psps.rs`) **have no caller left
  in `src/`**, which is the state E2 waits for;
- `every_census_in_the_cohorts_psps` (`census_cohort.rs`) has no non-test caller either —
  production calls the two halves, `each_census_in_the_cohorts_psps` and
  `the_censuses_as_one_cohort`. Keep it or delete it with a reason, but do not leave it undecided.

**E3 — the probes.** Both are already in `examples/`: `ng_census_read_vs_psp.rs` and
`ng_census_locus_spans.rs`. So E3 is *checking* rather than landing: that they build, that E2's
deletions do not break them, and that the first still writes its own census to a scratch file of
its own. **`examples/ng_census_locus_spans.rs` is one of the four files `cargo fmt` would rewrite
in the baseline** — leave that alone unless E4 decides to tidy it, and if it does, say so.

**E4 — words.** The subcommand docs for the three commands and the repair are current as of D1.
What is not: `PROJECT_STATUS.md`'s pipeline line (line 66, `generate-psps → generate-census →
estimate-parameters → call-from-psps`), and the plan's own §1 diagram is already right. The
implementation report for the whole plan is E4's too.

---

## How the work is run — the loop, unchanged

Per step: **implement → review → apply fixes → mutate → commit**, one step at a time, following
`ai/skills/plan-driven-implementation/SKILL.md`; every chat reply follows
`ai/skills/reporting-in-chat/SKILL.md`.

**The review arrangement.** One read-only agent per step over grouped categories, **forbidden to
edit any file, to write scratch files, or to run `cargo`**, asked to *name* mutations as exact
old→new text. The orchestrator runs them one at a time from a backup, proving each restore with
`diff`. Briefs from this milestone are in `tmp/b/d1_review_agent.md`, `d2_review_agent.md` and
`d3d4_review_agent.md`; the mutation harnesses are `tmp/b/d{1,2,3d4}_mutate.py` with drivers beside
them.

**Ten traps this plan has hit, all avoidable:**

- **A test that compares a file with itself proves nothing.** D1's whole suite passed with the write
  removed; D2's skip rule then made seven more tests compare the walk's bytes with themselves. The
  idiom that fixes it is in `regenerate_census/tests.rs`: `empty_every_trailer`, and
  `corrupt_a_block_of` for a psp that must fail on its records.
- **A mutation that does not compile, or that changes nothing, is indistinguishable from one nothing
  catches.** Two of D2's thirteen were faulty — one did not compile, one moved a block back where it
  already was — and the only reason it was visible is that the driver prints the apply step and the
  compiler's output beside each test result.
- **A background run that never started looks exactly like one that just started.** A driver of mine
  died on a syntax error and I reported it as running. Check the log grew, the first step reported,
  and the process is alive.
- **Never edit the tree while a gate or a mutation run is building it**, and **never run mutations
  while a review agent is reading the worktree**.
- **`tmp/b/gate.sh` is not executable**: run it as `nohup bash tmp/b/gate.sh <label> >
  tmp/b/<label>_gate.txt 2>&1 &`, and poll for `== fmt: files ==`.
- **`cargo fmt` rewrites four files that were already unformatted before this branch.** Run it, then
  `git checkout -- examples/ng_call_cohort_end_to_end.rs examples/ng_call_from_psps_cost.rs
  examples/ng_census_locus_spans.rs src/ng/psp/block.rs`.
- **The container runs as root**, so a file made read-only is still writable inside it.
- **Nothing counts a psp's block bytes.** `trailer_bytes_read` counts trailer bytes only, so a claim
  of the form *these records are not read* cannot be pinned; what can be pinned is that no rebuild
  pass happened, by corrupting the blocks. See the open item below.
- **Stale dev-container VMs starve the build** — each claims 16 GB; `container ls` shows them and
  `container stop <id>` clears them (the owner has approved stopping them).

---

## The gate, as a set

`bash tmp/b/gate.sh <label>` runs the four gates and prints each as a list to compare, never a
count. It takes 8-12 minutes.

| gate | at `0f4b0646` |
|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | **6,738 lib tests pass**, 0 failed, 15 ignored; 20 of 21 targets green, the one failure the baseline's `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | **11 errors of 5 kinds in 6 files** — `bam/alignment_input.rs`, `window_coverage/accumulator.rs` ×5, `cohort_merge/build.rs` ×2, `cohort_merge/serial.rs`, `psp_writer_line.rs`, `window_coverage/production_parity.rs` |
| `cargo check --all-targets --keep-going` | **4 examples** do not compile: `ng_call_cohort_end_to_end`, `ng_candidate_selection_probe`, `ng_cohort_merge_parallel_cost`, `ng_cohort_merge_real_cost` |
| `cargo fmt --check` | the **4 files** above |

**None of the four is this branch's**, and every step is judged against this set rather than against
green.

**The byte-identity oracle.** The parameters file fitted over `a_walked_cohort` in
`estimate_parameters/tests.rs` has been **20,741 bytes** and byte-identical since before C1
(`tmp/c1_before.toml`). A temporary test that writes `file.to_toml()` to a scratch path, run and
then deleted, re-checks it. Anything that changes it is a finding.

---

## Open, and the owner has not ruled on these

- **A block-bytes counter.** Nothing counts what a psp's blocks cost to read, so *a skipped psp's
  records are not read* cannot be tested — a mutation that reads them and discards the result passes
  everything. A counter beside `trailer_bytes_read` is a few lines. **My recommendation: add it in
  E if E touches the reader, and leave it otherwise.**
- **A three-sample psp fixture.** With two, a run that pressed on past a failure and rebuilt *later*
  samples is invisible, and the error doc's sixty-sample claim rests on the loop returning on the
  first failure. **Recommendation: leave it.**
- **A psp that will not open stops the cohort at the first file**, so spec §4.1's *every sample is
  examined* holds for stale censuses and not for unopenable files. **Recommendation: leave it.**
- **`cargo doc` is red**: 40 unresolved intra-doc links, none from this branch, and not in the gate.
  A separate cleanup, then add `cargo doc` to the gate.
- **`OpenPspCohort` is still in `psp_caller.rs`**, with eight users outside it. A pure move.
- **The seven-step preamble is duplicated** between `estimate-parameters` and `regenerate-census`,
  about 55 lines including two argument-identical catalog-check wrappers. Both are pinned by
  separate oracles, so a divergence would be caught rather than silent.

---

## Building, and where scratch goes

macOS with Apple's `container`: **every `cargo`, `python3` and `sed -i` goes through
`/Users/jose/devel/pop_var_caller-census-vs-psp/scripts/dev.sh`**. The shell's working directory
does not persist between commands and a new session's default directory may be the main checkout —
use absolute paths. Scratch goes under the worktree's `tmp/`, never the harness's own scratchpad.

Commit as `<type>(ng): E1 — <title>`, ending with
`Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`. Flip the step's `☐` to `✅` in the plan in
the same commit, with the code, the tests and both report sections. **Pause at Checkpoint E**, which
is the end of the plan — and at the end, per the plan-driven skill, also verify the build natively
on the host.
