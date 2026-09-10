# Handoff — the census lives inside the psp: D1 is committed, D2 to D4 remain

**Written 2026-09-10**, on branch `census-vs-psp-perf`, in the worktree
`/Users/jose/devel/pop_var_caller-census-vs-psp`. **The tree is clean at `34359259`.** This exists
so a fresh session can pick the plan up without re-deriving anything.

---

## Where the plan stands

**Milestone C is complete and D1 is committed.** The plan is
[`psp_census_pair.md`](../../ng/impl_plan/psp_census_pair.md); the design is spec
[`psp_census_pair.md`](../../ng/spec/psp_census_pair.md).

| commit | what |
|---|---|
| `4522a3e0` | C5 — a census recorded under other settings is refused, naming the setting |
| `93fb1ac3` | the spec's §4.2 row and the plan's C5 text say what the fit compares (owner ruled: *"Yes, edit the documents"*) |
| `34359259` | D1 — `generate-census` becomes `regenerate-census`, which replaces each psp's trailer |

**Two things are owed on D1 before D2 starts**, and neither is a code change so far as anyone
knows:

1. **The step's review was never applied.** One read-only agent was launched over D1's diff and the
   session ended before it reported, so its findings are lost. The brief it was given is
   [`tmp/b/d1_review_agent.md`](../../../../tmp/b/d1_review_agent.md) — re-run one agent with it,
   against `git show 34359259` rather than the working tree, and fix forward in its own commit.
2. **D1's mutations were not run.** The harness is `tmp/b/c5_mutate.py` plus
   `tmp/b/c5_mutations.sh`; copy them for D1's own list, which the review is asked to name.

---

## What D2, D3 and D4 need, beyond what the plan says

Read the plan's own text for each step first. What follows is what this session learned that the
plan does not say.

**D2 — fresh psps skipped.** Spec §8's rule is freshness *whole*: all three of §4.2's causes, since
this command rebuilds the selection anyway. Every piece now exists and none of it did when the plan
was written:

- the head's two causes: `what_the_heads_say_about_every_census_in_a_cohort` (`census_freshness.rs`);
- the third: `CensusPlan::recording_terms` (`gatherer.rs`) against each census's own recorded
  settings, through `CensusVerdict::of_recorded_settings` or
  `what_the_run_says_about_every_census_in_a_cohort`;
- each census read cheaply, without its records: `each_census_in_the_cohorts_psps`
  (`census_cohort.rs`) decodes a census's header and directory and nothing else.

**The order D1 leaves behind is what D2 has to fit into**: the cohort is opened, the reference and
the catalog are checked against the psps, the selection is planned, and *then* the cohort is closed
and the psps are rewritten one at a time. A skip decision needs the plan, so it belongs after the
plan is built — and the psps are still open at that point, which is where the head's read wants
them. `SampleCensusOutcome` and `CensusReport::lines` will need a skipped count; the report's first
line says *regenerated N censuses* today.

**D3 — a stopped run costs only what is left.** It is D2's property, tested: regenerate three, make
the third psp unreadable after two succeed, then readable again, and the second run does the third
alone. Nothing in D1 stands in the way — each sample is one `census_from_psp` and one
`replace_trailer`, and a failure returns without touching the rest.

**D4 — the whole-file oracle.** D1 already asserts the *trailer* the walk wrote equals the trailer
the rebuild writes, on two fixtures, plus that the header and the record count are unchanged. What
D4 adds is the whole file, `cmp`-identical, on a copy — and the same on real reads through
`scripts/ng_fit_stage_end_to_end.sh`, which is step E1's own work and still calls `generate-census`.

---

## How the work is run — the loop, unchanged

Per step: **implement → review → apply fixes → commit**, one step at a time, following
`ai/skills/plan-driven-implementation/SKILL.md`; every chat reply follows
`ai/skills/reporting-in-chat/SKILL.md`.

**The review arrangement.** One read-only agent over grouped categories (reliability and errors,
idiomatic and smells, naming and module structure, prose against code, test strength), **forbidden
to edit any file, to write scratch files, or to run `cargo`**, and asked to *name* mutations as
exact old→new text. The orchestrator runs them, one at a time, from a backup, proving each restore
with `diff`.

**Six traps this plan has hit, all avoidable:**

- **Never edit the tree while a gate or a mutation run is building it.** Twice in this session a
  doc-comment edit landed mid-gate and the gate had to be discarded and re-run.
- **Never run mutations while a review agent is reading the worktree.**
- **A mutation that does not apply looks exactly like one nothing catches.** The harness asserts its
  match count; never silence the apply step.
- **`tmp/b/gate.sh` is not executable.** `nohup tmp/b/gate.sh d1 &` fails with *permission denied*,
  the launcher exits 0, and the only sign is one line in the log. Run it as
  `nohup bash tmp/b/gate.sh <label> > tmp/b/<label>_gate.txt 2>&1 &`.
- **`cargo fmt` rewrites four files that were already unformatted before this branch.** Run it, then
  `git checkout -- examples/ng_call_cohort_end_to_end.rs examples/ng_call_from_psps_cost.rs
  examples/ng_census_locus_spans.rs src/ng/psp/block.rs`. None of them holds this branch's work.
- **The container runs as root**, so a file made read-only is still writable inside it. A test that
  wants a write to fail cannot get there that way; D1 pins the message of the refusal instead.

**Stale dev-container VMs starve the build** — each `./scripts/dev.sh` VM claims 16 GB; `container
ls` shows them and `container stop <id>` clears them (the owner has approved stopping them).

---

## The gate, as a set

`bash tmp/b/gate.sh <label>` runs the four gates and prints each as a list to compare, never a
count. It takes 8-12 minutes; launch it with `nohup … &` and poll for `== fmt: files ==`.

| gate | at `34359259` |
|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | **6,731 lib tests pass**, 0 failed, 15 ignored; 20 of 21 targets green, the one failure the baseline's `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` |
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

## Still open, and the owner has not ruled on these

Recommendation for all three: leave them.

- **A psp that will not open stops the cohort at the first file**, so spec §4.1's *every sample is
  examined* holds for stale censuses and not for unopenable files.
- **`cargo doc` is red**: 40 unresolved intra-doc links, none from this branch, and not in the gate.
  A separate cleanup, then add `cargo doc` to the gate.
- **`OpenPspCohort` is still in `psp_caller.rs`**, with seven users outside it. A pure move, worth
  doing early in whatever step next touches that file.

---

## Building

macOS with Apple's `container`: **every `cargo`, `python3` and `sed -i` goes through
`/Users/jose/devel/pop_var_caller-census-vs-psp/scripts/dev.sh`**. The shell's working directory does
not persist between commands and a new session's default directory may be the main checkout — use
absolute paths. Scratch goes under the worktree's `tmp/`, never the harness's own scratchpad.

---

## Reports and commit convention

Milestone D's report is
[`ng_psp_census_pair_milestone_d_2026-09-10.md`](ng_psp_census_pair_milestone_d_2026-09-10.md); its
D1 section is the context for D2. **The milestone has no review report yet** — D1's review never
landed — so the next session creates
`doc/devel/reports/reviews/psp_census_pair_milestone_d_2026-09-10.md`.

Commit as `feat(ng): D2 — <title>`, ending with
`Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`. Flip the step's `☐` to `✅` in the plan in
the same commit, with the code, the tests and both report sections. **Pause at Checkpoint D** — the
plan marks it, and the owner reviews there.
