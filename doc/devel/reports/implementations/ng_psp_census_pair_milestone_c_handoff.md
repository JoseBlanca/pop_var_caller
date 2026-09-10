# Handoff — the census lives inside the psp, step C5, then hand back

**Written 2026-09-10**, on branch `census-vs-psp-perf`, in the worktree
`/Users/jose/devel/pop_var_caller-census-vs-psp`. **Milestone C (C1–C4) is committed and the tree is
clean.** This exists so a fresh session can pick the plan up without re-deriving anything.

---

## What is left in the next session's scope

**One step, C5, and then a hand-back before Milestone D.** It was added to the plan at Checkpoint C
by the owner's ruling of 2026-09-10 — read its full text in
[`psp_census_pair.md`](../../ng/impl_plan/psp_census_pair.md) (the step after C4). In one sentence:
**when a census was built under settings other than this run's, refuse to fit, and say which
setting differs and what fixes it.**

The owner's words: *"that's what we should do, stop and report the problem."*

---

## What C5 has to get right, already established

1. **Why the setting has to be named.** A census records seven settings (`SELECTION_FIELDS`,
   `src/ng/parameter_estimation/joint/census.rs`): selection seed, reference digest, analysed region
   set, repeat catalog build settings, STR routing criteria, generic target position count, STR
   per-stratum cap. A mismatch in different ones has **opposite fixes**: a different reference or
   catalog means the user passed the wrong file and should rerun with the right one, which rebuilds
   nothing; a different seed, budget or cap means this build picks positions differently and only
   regenerating fixes it — a quarter of an hour a sample at whole-genome scale. A message that can't
   say which leaves the user to guess, and a wrong guess costs hours and ends refused again.
2. **What happens today, measured.** The fit compares only a digest of the *kept positions*
   (`fit_a_cohort`, `src/ng/run/census_fit.rs`). At C4, a cohort whose censuses were rebuilt under
   **half the shipped position budget** fitted without a word on the fixture, because its 600-base
   contig holds fewer ordinary positions than either budget, so both selections kept the same set.
   That cohort is C5's first test: it must now be refused, naming *generic target position count*.
   The fixture recipe is in C4's test
   `a_cohort_written_against_another_selection_is_told_both_ways_out`
   (`src/pop_var_caller_exp/estimate_parameters/tests.rs`): `a_census_plan_over_selecting`, then
   `census_from_psp`, `write_census(.., None, ..)` and `replace_trailer` per psp.
3. **The pieces are in hand.** `CensusPlan.terms: SelectionTerms` (`src/ng/run/gatherer.rs`) is the
   run's; `CohortCensusEvidence::terms()` gives a census's `RecordingTerms`, whose `selection` field
   is a `SelectionTermsDigest`; `SelectionTermsDigest::of(&plan.terms).first_disagreement(..)`
   returns the field name. The kept-positions digest stays as a backstop.
4. **A wrong `--reference` should be caught before any selection is rebuilt.** `call-from-psps`
   already compares the reference with every psp header (`refuse_a_file_against_another_reference`,
   private in `src/ng/run/psp_caller.rs`); `estimate-parameters` never does. Spec §6 also claims
   `--catalog` is checked against the header's catalog — it is not; C5 should make that true.
5. **A second message names no action.** A cohort walked by two builds with different selection
   constants passes the freshness judgement and is refused while its censuses are read
   (`CohortCensusEvidence::new`) as *samples A and B disagree on selection seed; they did not record
   the same thing*. **Ordering trap:** that refusal happens *before* the reference is read, while the
   comparison against the run needs the reference read and the selection rebuilt. Decide deliberately
   where each lives; the plan text prefers the cross-sample case as rows of C3's report
   (`CensusesToRegenerate`, `src/ng/run/census_freshness.rs`).
6. **The command name.** Both existing messages say `regenerate-census` through one constant,
   `THE_COMMAND_THAT_REBUILDS_A_CENSUS` in `census_freshness.rs`. The command does not exist until
   D1. At D1 the new subcommand's name is defined *from* this constant, never the reverse — `ng`
   imports nothing from `pop_var_caller_exp` outside tests. Plan step D1's text does not mention it.

---

## Still open from Checkpoint C — the owner has not ruled on these

Recommendation for both: leave them. Ask only if C5 makes one of them matter.

- **A psp that will not open stops the cohort at the first file**, so spec §4.1's *every sample is
  examined* holds for stale censuses and not for unopenable files.
- **`cargo doc` is red on this tree**: 40 unresolved intra-doc links, none from this branch, and not
  in the gate. Recommendation: separate cleanup, then add `cargo doc` to the gate.

And one suggestion from Checkpoint B not yet done: moving `OpenPspCohort` out of `psp_caller.rs`
into its own module (it now has five users outside it). A pure move; worth doing early in D.

---

## How the work is run — the loop, unchanged

Per step: **implement → review → apply fixes → commit**, one step at a time. The plan-driven skill is
`ai/skills/plan-driven-implementation/SKILL.md`; every chat reply follows
`ai/skills/reporting-in-chat/SKILL.md`.

**The review arrangement.** One read-only agent over grouped categories (reliability + errors,
idiomatic + smells, naming + module structure, refactor safety), forbidden to edit files or run
`cargo`, asked to *name* mutations as exact old→new text. The orchestrator runs them, one at a time,
from a backup, proving the restore with `diff`. **Every step of Milestone C had at least one test
that looked sound and could not fail**; the mutations are where the value is.

**Four traps this milestone hit, all avoidable:**

- **Never run mutations while a review agent is reading the worktree.** C3's reviewer read a file
  with a mutation applied. Finish mutations first, or give the agent a copy.
- **A mutation that does not apply looks exactly like one nothing catches.** The harness asserts its
  match count, but sending the apply step to `/dev/null` in the same command hid the assertion once.
  Never silence the apply step.
- **Never `git checkout --` a file holding uncommitted work.** Restore mutations from their backup.
- **Stale dev-container VMs starve the build.** Each `./scripts/dev.sh` VM claims 16 GB. On
  2026-09-10 five abandoned ones (Aug 28 – Sep 9) left 650 MB free and the gate stalled for 20
  minutes; `container ls` shows them and `container stop <id>` clears them (the owner approved
  stopping them). One rustc internal compiler error during the gate was transient and passed on
  re-run.

---

## The gate, as a set

`tmp/b/gate.sh <label>` runs the four gates and prints each as a list to compare, never a count. It
takes 8–12 minutes (four container builds); launch it with `nohup … &` and poll the summary file.

| gate | at `9733f6cb` (end of Milestone C) |
|---|---|
| `cargo test --lib --bins --tests --all-features --no-fail-fast` | **6,723 lib tests pass**, 15 ignored; 20 of 21 targets green; the one failure is `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | **11 errors of 5 kinds in 6 files** (`bam/alignment_input.rs`, `window_coverage/accumulator.rs` ×5, `cohort_merge/build.rs` ×2, `cohort_merge/serial.rs`, `psp_writer_line.rs`, `window_coverage/production_parity.rs`) |
| `cargo check --all-targets --keep-going` | **4 examples** do not compile |
| `cargo fmt --check` | **4 files** |

**`cargo fmt` reformats those four files too**: run it, then
`git checkout -- examples/ng_call_cohort_end_to_end.rs examples/ng_call_from_psps_cost.rs
examples/ng_census_locus_spans.rs src/ng/psp/block.rs` before committing.

**The byte-identity oracle.** The parameters file fitted over `a_walked_cohort` in
`estimate_parameters/tests.rs` has been 20,741 bytes and byte-identical since before C1
(`tmp/c1_before.toml`). A temporary test that writes `file.to_toml()` to a path from an environment
variable, run and then deleted, re-checks it. C5 should leave it unchanged for a fresh cohort.

---

## Building

macOS with Apple's `container`: **every `cargo`, `python3` and `sed -i` goes through
`/Users/jose/devel/pop_var_caller-census-vs-psp/scripts/dev.sh`**. The shell's working directory does
not persist between commands, and a new session's default directory is the main checkout — use
absolute paths. Scratch goes under the worktree's `tmp/`.

---

## Reports and commit convention

Append C5 to both Milestone C reports:
[`ng_psp_census_pair_milestone_c_2026-09-10.md`](ng_psp_census_pair_milestone_c_2026-09-10.md) (its
Checkpoint C section is the context for C5) and
[`../reviews/psp_census_pair_milestone_c_2026-09-10.md`](../reviews/psp_census_pair_milestone_c_2026-09-10.md).

Commit as `feat(ng): C5 — <title>`, ending with
`Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`. Flip C5's `☐` to `✅` in the
plan in the same commit, with the code, the tests and both report sections.
