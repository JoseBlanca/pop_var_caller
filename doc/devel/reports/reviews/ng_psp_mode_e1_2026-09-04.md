# Code Review: ng_psp_mode_e1
**Date:** 2026-09-04
**Reviewer:** rust-code-review skill (orchestrator; two sub-agents in isolated worktrees, seven category checklists between them)
**Scope:** step E1's uncommitted diff — opening a cohort of psps and refusing one that cannot be called, on `d2c7113e`
**Status:** Request-changes (all applied — see the fix report)

---

### 1. Scope

- Reviewed: the working-tree diff of plan step E1, exported as `tmp/review_2026-09-04_e1/e1.patch`.
- In-scope: [psp_caller.rs](../../../src/ng/run/psp_caller.rs) (new), and the changes to [run/mod.rs](../../../src/ng/run/mod.rs), [callers.rs](../../../src/ng/run/callers.rs) and [read_groups.rs](../../../src/ng/read/input/read_groups.rs).
- Categories (7): reliability, errors, defaults, naming, module_structure, idiomatic, smells.

### 2. Verdict

**Request-changes**, on more than the usual. The two agents between them found **10 Major**, and the most valuable of them were not defects in what was written but **checks that were missing**: the file-descriptor refusal spec §7.1a asks for by name, direct mode's catalog-against-reference refusal, and the psp header's own whole-assembly digest — the field's one documented consumer, which nothing was consuming.

One finding is a **spec-against-spec conflict** and is carried to the checkpoint rather than settled here: §6.2 asks for a refusal that the psp format's own validator declares legal.

### 3. Execution status

- Both agents detached at `d2c7113e`, applied the patch, verified the marker file, and restored their trees.
- Agent mutations: 4 + 5 run, 1 + 5 survived.
- Orchestrator: 10 mutations before the review, all killed; 8 after it, 5 killed and 3 survived; the two real survivors closed and re-killed. **18 in total, 17 killed, 1 changed no reachable behaviour.**

### 4. Open questions — both carried to Checkpoint E

1. **§6.2 asks for a refusal the format declares legal.** Two entries of one sample sharing an `@RG ID`: the spec's reason is that such a table cannot be renumbered without guessing, and this format guesses nothing — identity is the walk-local number, which is the entry's position. `psp/header.rs`'s validator says a sample sequenced across files may legally carry two such entries, and direct mode calls that cohort. **Not refused, with the reasoning recorded at the refusal; the spec clause needs amending.**
2. **§6.2's by-name parameters match is a count.** `RunParameters` carries no names, so the by-name half cannot live here; it belongs at F1, where the parameters file meets the cohort's sample list. **Recorded at the call site; the plan's E1 entry needs amending.**

### 5. Top 3 priorities

1. **Three refusals missing**, one of them named in the spec's own §7.1a.
2. **Every per-file refusal pinned only on a one-file cohort** — a mutant that checked `psps[0]` and stopped passed all 490 tests.
3. **The within-file read-group order untested** — reversing it passed all 490 tests, and it silently puts every observation on the wrong calibration.

### 6. Findings

#### Major

- **No file-descriptor refusal.** Spec §7.1a: *"a check at construction, beside the header checks of §6.2"*. Direct mode has one; psp mode opened every file first and would die at file 249 on a macOS default limit, blaming an innocent path.
- **Direct mode's catalog-against-reference refusal has no psp-mode counterpart.** A catalog built on another build of the assembly puts every repeat tract at the wrong position, and every segment the run loops over is drawn from it — refused in one mode and not the other, over the same inputs.
- **The psp header's whole-assembly digest was never read**, though its own doc gives this check as its purpose.
- **Every per-file refusal was pinned on a one-file cohort.** Mutation: the run loop `.take(1)` — 490 tests green.
- **The within-file read-group order was untested.** Mutation: walk each file's table backwards — 490 tests green, and `remap[0]` becomes the run number of the file's *last* group.
- **`of_merged_tables`' guard did not check what its own `# Panics` promises.** It checked coverage and never *whose*: a table where one sample's entry claims another's read group passed all three assertions — verbatim the symptom the panic message names.
- **Library and experiment both came back `Declared`**, including names a walk invented. Direct mode always marks the experiment `Synthesized`, with the reason written beside it.
- **The duplicate-`@RG ID` refusal fires on a legal table** — open question 1.
- **The by-name parameters match is a count** — open question 2.

#### Minor

`sample_count` was inserted into the middle of `into_sources`' doc comment, so rustdoc gave each the other's text; a doc claim that the caller's segmentation cannot be checked, when it is checked three lines below; "the header is the second thing a reader reaches" (it is the third); "cannot be called from here" of the ground assembly, which would in fact compile — the constraint is dependency direction; "§12.3's byte-identity oracle" (§12.3 is mode equivalence); "123 kB a cursor" quoted for a stage that holds no cursor; a duplicated per-sample remap vector, one field away from the hazard `of_merged_tables` panics about; `of_merged_tables` untested in its own module, with its whole guard removable and everything still green; and `#[allow(dead_code)]` where this crate already uses `#[expect]` so the first real reader turns the line into a compile error.

### 7. What holds up

- **Every refusal names what a person can act on** — the sample, the file, the field, the contig, or the two numbers — and each has a test that reads it back.
- **The two-move opening is the right shape**, and the second agent proved the layering constraint behind it rather than taking the doc's word: reaching down for the command layer compiles, so the argument is about dependency direction and now says so.
- **The order of the refusals is defensible and in one place better than direct mode's**: psp mode needs nothing from its files to compare its two views of the reference, so that check is hoisted above the opens where direct mode must defer it.

### 8. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib                  # 6,150 passed
./scripts/dev.sh cargo test --lib 'ng::run'        # 496 passed
bash tmp/e1_mutations/run4.sh; bash tmp/e1_mutations/run5.sh
```
