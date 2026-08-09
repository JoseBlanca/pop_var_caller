---
name: rust-code-review
description: Use this skill whenever the user asks for a code review, audit, critique, or quality check of Rust code (.rs files, Cargo crates, snippets, PRs, or diffs). Trigger on phrases like "review my Rust code", "audit this crate", "is this idiomatic Rust", or when Rust source is shared with any request implying quality feedback.
---

# Rust Code Review

You are performing a professional, uncompromising code review of Rust code. Quality is the highest priority — be precise, specific, and direct. Vague praise is forbidden; every comment must point to a concrete location and propose a concrete change.

The review is split across focused per-category checklists in `ai/skills/rust-code-review/code_review/`. **You are the orchestrator**: you triage which categories apply to the scope, dispatch one sub-agent per category in parallel — **each in its own git worktree** (step 6) — then synthesize their findings into a single report. Per-category rules are not duplicated in this file — read each category file when dispatching the corresponding sub-agent.

## Review principles (must always hold)

- **Correctness over style.** Prioritize behavioral correctness and failure transparency over formatting or personal preferences.
- **No silent assumptions.** When the code leaves something unspecified — invariants on inputs, call-site guarantees, threading context, whether a collection is sorted, whether a value can be zero or empty — do not silently pick an answer and review as if it were fact. State the assumption explicitly in the finding and lower its severity until it can be verified.
- **Evidence-first findings.** Do not invent file paths, line numbers, logs, or command results. If you cannot verify, label as "Needs verification".
- **Actionability.** Every non-trivial finding must include a concrete fix (diff/snippet, test, or refactor step), or — if the right fix depends on intent the reviewer cannot infer — a specific question whose answer would determine the fix.
- **Scope discipline.** Review what was asked. For a diff or PR, focus on changed lines and their direct callers/callees; flag pre-existing issues in untouched code separately under "Out of scope observations". Exception: pre-existing Blocker-severity issues (security, data loss, undefined behavior) are raised under Findings regardless of scope.

The severity rubric and per-finding format are defined in `ai/skills/rust-code-review/code_review/_finding_format.md`. Read it once at the start of every review — both you (for synthesis) and every sub-agent you dispatch will follow it.

## Review procedure

You run each numbered step once per review. Only step 6 fans out into parallel sub-agents.

### 0. Read `PROJECT_STATUS.md`

Read `PROJECT_STATUS.md` at the project root to orient on what the project is, the current focus, and the in-scope feature's prior artefacts (plan, implementation report, prior reviews, prior fix-applied reports). See *Project status protocol* below for what to read and what to update at the end.

### 1. Establish scope

Determine whether the review covers a full crate, a diff/PR, or a snippet. State the scope. For diffs, identify changed files and the direct callers/callees of changed items; these define the in-scope surface.

### 2. Inventory

List files, public API surface, error types, concurrency primitives, and external dependencies. Detect category triggers: does this code use `unsafe` / `async` / `Arc` / `Mutex` / atomics / channels? Does it have a public API? Is there a `Cargo.toml`? Is it a parser / validator / security boundary? Is it on a hot path?

### 3. Run verification commands

If a real execution environment is available, run, in the project's container per `CLAUDE.md`:

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo doc --no-deps`
- `cargo audit`

Quote actual output verbatim. **Never simulate, paraphrase, or guess at command output.** If commands cannot be run, list them under "Commands the author must run" and proceed with static review only. Any real failure is at least **Major**; correctness-impacting failures are **Blocker**. The verbatim output is passed into each sub-agent's prompt at step 6 so they do not re-run.

### 4. Determine intent

Before judging whether the code works, establish what it is meant to do — its purpose, its contract with callers, the inputs it is built to handle, the invariants it maintains. Correctness is meaningless without intent. If intent is unclear from names, types, tests, and docs, file it as a finding against documentation or naming at Minor severity (Major if the ambiguity could plausibly cause misuse).

Write a one-paragraph "domain intent" summary; it is passed into each sub-agent's prompt.

### 5. Triage categories

Decide which per-category checklists apply. Each lives at `ai/skills/rust-code-review/code_review/<category>.md`.

| Category | Apply when |
|---|---|
| `reliability` | Always. Snippets without test files: still flag missing tests. |
| `errors` | Always. |
| `naming` | Always. |
| `defaults` | Scope contains public API, configuration, or any default-acting value. Skip pure-internal snippets with no parameters. |
| `idiomatic` | Always. |
| `refactor_safety` | Always. |
| `module_structure` | Scope spans multiple files (crate / PR / module-level review). Skip pure single-file snippets. |
| `unsafe_concurrency` | Code uses `unsafe`, `Arc`, `Mutex`/`RwLock`, atomics, channels, `async`, or thread spawning. Skip otherwise. |
| `smells` | Always. |
| `tooling` | Scope is a crate (has `Cargo.toml`). Skip pure snippets. |
| `extras` | Scope contains parsers/validators/security boundaries; accepts untrusted input; produces stable output; is on a hot path; is a public crate; or is a PR (for "Diff matches stated intent"). Apply per item. |

When in doubt, dispatch — a sub-agent that finds nothing applicable writes `No findings.` and is cheap.

### 6. Dispatch sub-agents in parallel — **each in its own worktree**

Create the scratch directory: `tmp/review_<YYYY-MM-DD>_<scope-slug>/` (append `_v<N>` if it already exists). The slug is a short kebab-case identifier of the reviewed module or PR (e.g. `gvcf_parser`, `pr-142`).

Project rule: scratch space is project-local `tmp/`, never `/tmp`. Add `tmp/` to `.gitignore` if it is not already covered by the existing target ignores.

For each selected category, dispatch a `general-purpose` sub-agent **in parallel** — issue a single message with multiple Agent tool calls — and pass **`isolation: "worktree"`** on every one.

#### Why the isolation is mandatory, not an optimization

A reviewer that only reads cannot tell a load-bearing line from a decorative one. The findings that matter come from *changing* the code and watching what fails — and the moment two agents do that in one checkout they are editing the same files at the same time. Measured, on one 9-agent review of a shared tree:

- agents reported their edits overwritten mid-experiment, `src/` reverted under them three times, and one build reading a half-written file;
- one agent's baseline failed on **another agent's** marker, which it then had to diagnose;
- five of the nine detected the interference and retreated into private worktrees of their own — so more than half the fan-out paid the isolation cost anyway, late, after wasted work;
- every mutation result was untrustworthy in principle (a "green" run can be another agent's revert landing first, a "red" one their mutation), so the orchestrator had to re-verify the decisive findings **serially** before accepting them.

The next milestone's review, over comparable code with 5 agents and a worktree each, had zero collisions, every result first-hand and nothing needing re-verification. Isolation is cheaper than the bookkeeping its absence forces.

#### The worktree starts on `main`, not on your branch — every prompt needs a step 0

**`isolation: "worktree"` hands each agent its own tree and checks it out at `main`.** Probed
during Milestone D's review: the agent's `HEAD` came back exactly `git rev-parse main`, and
`src/ng/locus_generation/pileup/` held 11 files there against 13 on the branch — the two missing
ones being the generator itself. **The code under review was absent from every agent's tree.**

The failure is silent and expensive. An agent handed its scope by absolute path finds those files
in *your* checkout, because that is the only copy on the machine where they exist — so a fan-out
meant to isolate five mutation-heavy agents converges all five on the author's working tree. That
happened once, and the whole review had to be re-run.

So **every agent prompt begins with a step 0**, before anything else:

> 0. **Re-point your worktree at the commit under review.** Your tree was created from `main` and
>    does not contain this branch's code. Run `git checkout --detach <sha>` in your own worktree
>    root, then `ls` <two files that exist only on the branch>. If either is missing, **stop and
>    report** — do not review from another checkout.

Verified on a worktree created off `main` exactly as the harness makes them: the branch-only files
were absent before and present after, the tree stayed clean, and no fetch was needed since the
worktrees share one object store. A second mechanism stacks on top: an **unchanged** worktree is
auto-cleaned, so an agent resumed after that lands in a real checkout — one more reason the check
is a hard stop rather than a warning.

#### What isolation changes in the prompt

Three things, and all three must be said explicitly or the agent will get them wrong:

1. **Where it builds — its own copy of the build script, by that copy's absolute path.** The agent
   is *in* its worktree, and this project's `scripts/dev.sh` derives `PROJECT_DIR` from the
   script's own location: **the agent's copy mounts and builds the agent's tree**, while
   `<main-checkout>/scripts/dev.sh` builds the main checkout and ignores the caller's directory.
   Instruct the agent by its own worktree's absolute path. (Measured: an agent's own `dev.sh`
   compiles that worktree in ~31 s and leaves ~1.3 GB in its `target-container`. An earlier note
   claiming `dev.sh` "cannot build a worktree" was describing an invocation of the *main* copy —
   it is wrong, and it is what misdirected one fan-out.)
2. **Where its findings go.** Its worktree is temporary and is cleaned up. The findings file must be written to an **absolute path in the main checkout's** `tmp/review_<date>_<slug>/`, which is outside the worktree and survives. Each agent writes its own file, so there is no contention.
3. **What it may now do.** Say plainly that it has its own tree and should mutation-test aggressively — the isolation only pays if the agents use it. An agent that still reviews by reading has cost you a worktree for nothing.

**Clear the worktrees afterwards, and lift the evidence out first.** A dirty worktree is evidence
only until it is deleted, so extract any probe or candidate fix you want to keep before pruning:
`git worktree remove --force <path>` then `git branch -D worktree-agent-<id>` (the branch pointer
survives the removal).

Do **not** ask isolated agents to leave the tree clean or to revert their experiments. That instruction is for a shared checkout; here it wastes their effort on a tree that is about to be discarded, and a dirty worktree is evidence you can still inspect.

Each sub-agent prompt:

> Run the **<category>** checklist on the following Rust code review scope.
>
> **Scope:** <full crate / PR diff / snippet>
> **Domain intent:** <one paragraph from step 4>
> **In-scope files (full paths):** <list>
> **Out of scope:** <list, with reasons>
> **Verification command output:** <verbatim quotes from step 3, or "not run, because …">
>
> **Instructions:**
> 0. **Re-point your worktree at the commit under review.** It was created from `main` and does
>    not contain this branch's code. Run `git checkout --detach <sha>` in your own worktree root,
>    then confirm <two branch-only files> exist. If either is missing, **stop and report** rather
>    than reviewing from another checkout.
> 1. Read `ai/skills/rust-code-review/code_review/<category>.md` for the rules to apply.
> 2. Read `ai/skills/rust-code-review/code_review/_finding_format.md` for the severity rubric and finding format.
> 3. Read each in-scope file.
> 4. Apply each rule and produce findings in the specified format.
> 5. **You are in your own isolated git worktree.** Build, mutate and revert freely *there* — nothing you do can collide with another agent, so mutation-test aggressively rather than reviewing by reading. Run the project's build tooling **from your own worktree's copy, by its absolute path** (`<your-worktree>/scripts/dev.sh …`): the script builds whichever tree it lives in, so the main checkout's copy would build the wrong code.
> 5a. **Before recording a mutation as a survivor, prove it changed behaviour.** A mutant that takes the same path on every fixture is not a finding, and reporting one sends the author to defend a hazard that does not exist. Report three numbers, not two: mutations run, survived, and changed-no-behaviour.
> 5b. **Also challenge the tests that already exist.** For each, name what in its fixture makes the asserted failure reachable; if nothing does, that is the finding. Ask of every assertion which *wrong* implementations also satisfy it — if several do, the fixture is not discriminating.
> 5c. **Tell me what I already know.** If the dispatch listed mutations already run and their outcomes, do not spend a round re-finding them; verify them if you doubt them, and say so.
> 6. Write findings to the **absolute path** `<main-checkout>/tmp/review_<date>_<slug>/<category>.md` — outside your worktree, which is temporary and will be discarded. If no findings apply, write only the line `No findings.`
> 7. Do not invent file paths, line numbers, command output, or behavior. Cite only locations you have read.
> 8. Stay within the category. Issues that belong elsewhere go under a `## Cross-category observations` heading at the bottom of your file.

Substitute `<category>` and the scope fields for each dispatch. Do **not** assign severity codes (B1, M1, Mi1, …) inside sub-agents — that happens at synthesis.

### 7. Collect findings and route cross-category observations

For each `tmp/review_<date>_<slug>/<category>.md`, in turn:

1. **Read the file.** If it contains only `No findings.`, mark the category as clean and continue.
2. **Tally findings** by severity. Preserve each sub-agent's per-finding text — you will paste it during synthesis (step 9), with a severity code prepended.
3. **Read the `## Cross-category observations` section at the bottom.** Sub-agents are instructed to note issues that belong to other categories there rather than filing them out-of-category. For each note, decide:
   - **Promote to a finding** if the issue clearly matches a rule in the destination category and was not already raised there. Add it to that category's tally; during synthesis, file it under the destination category's severity and cite that multiple categories surfaced it ("convergent finding").
   - **Merge** if the same issue is raised in multiple places (its own category's findings *and* one or more cross-category notes from other agents). File once during synthesis, citing every category that surfaced it. Convergent evidence raises confidence; do not duplicate the entry.
   - **Defer** if the note is genuinely out of scope (different module, different concern) — record it for "Out of scope observations" in the synthesized report.
4. **Sanity-check sub-agent output.** If a sub-agent appears to have skipped its scope, cited locations it did not read (a hallucination risk), or produced findings without concrete fixes, redispatch *that one category* with an explicit instruction to fix the gap. Do not redo the others.

The point of this step is to do the routing once, deterministically, before synthesis — not to leave it as an exercise during the writeup.

### 8. Verify the test challenge

The `reliability` sub-agent runs the "challenge tests" pass for every non-trivial function. Spot-check that it did so for the changed-code surface; if a non-trivial function was missed, supplement with the missing entries before synthesis.

### 8a. Verify the diff's own quantitative claims

**Every number in a changed doc comment, test comment or commit message that describes *this work's* reach is measured or it is wrong.** How many cases a fixture covers, how far apart two values are, what a mutation cost, how many tests were added — re-derive each one, and report it CHECKED-CORRECT or WRONG with the right value.

This is its own step and not a category because it is about the diff's prose rather than any one file, and because it is where the errors concentrate. Measured on one plan in this repo: **fifteen wrong numbers in a single milestone, every one the author's own claim about their own fixture, while about forty figures quoted from the design and research documents were all correct.** A figure copied from a source is usually right; a figure describing the author's own test is where to look.

Two habits make the step cheap:

- **Compute, do not read.** A claim that a test covers 2,875 cells is checked by running the fixture, not by re-reading the sentence. Several wrong numbers survived a report, a commit message and a chat summary because everyone re-read them.
- **Check the mechanism, not only the magnitude.** An explanation of *why* something behaves as it does is a claim, and a wrong one is worse than a wrong number: it sends the next reader hunting a symptom that does not occur. Reproduce the story. Two examples from one milestone: a comment attributing a fixture's insensitivity to one cause when restoring that cause makes the test fail on *correct* code, and a comment justifying a selection rule by a failure mode that the measurement shows scores 2,165 nats in the opposite direction.

### 9. Synthesize the unified report

Compose the report using the *Output format* below. Verdict, top 3, and "What's good" need the full picture and are produced by you. Assign severity codes during synthesis: `B1, B2, …` for Blocker; `M1, M2, …` for Major; `Mi1, Mi2, …` for Minor; Nits stay grouped without numbering.

Each finding is filed once. Issues that were merged in step 7 (raised by multiple categories) carry a `**Categories:** <a>, <b>` line citing every sub-agent that surfaced them — this is convergent evidence and should be visible to the author, not hidden behind deduplication.

Save to `doc/devel/reports/reviews/<module-slug>_<YYYY-MM-DD>.md` per the saving conventions below. Leave the per-category files in `tmp/` as an audit trail.

### 10. Update `PROJECT_STATUS.md`

After the review report is saved, update the in-scope feature's block in `PROJECT_STATUS.md` to point at it. See *Project status protocol* below for the rules.

## Output format

Produce the synthesized report in the following order. Use the section headings verbatim so the format is machine-readable.

### 1. Scope
- What was reviewed: full crate / PR diff / snippet.
- Reviewed against: commit hash, branch, or "as-provided".
- In-scope files (list).
- Deliberately out of scope (list, with reason).
- Categories dispatched (list, each with a one-line reason).

### 2. Verdict
Approve / Approve-with-changes / Request-changes.

### 3. Execution status
- Commands run, with exit code and one-line result.
- Commands not run, and why.
- Count of findings labeled "Needs verification".

### 4. Open questions and assumptions
Numbered. Each entry references the findings it affects. The author resolves these before responding to individual findings.

### 5. Top 3 priorities
The highest-impact fixes, with one-line rationale and pointer to the full finding.

### 6. Findings
Grouped by severity (**Blocker** → **Major** → **Minor** → **Nits**). Within each severity, ordered by confidence (High → Low), then by file. Nits collected into a single sub-section, not enumerated. Each finding follows the format defined in `ai/skills/rust-code-review/code_review/_finding_format.md`, with the severity code (B1, M1, …) prepended to the title — e.g. `B1: src/parser.rs:42 — Title`.

### 7. Out of scope observations
Pre-existing issues in untouched code, surfaced but not blocking. Each: file, brief description, suggested follow-up (separate PR or issue). Pre-existing Blocker-severity issues (security, data loss, UB) appear under Findings instead, marked "pre-existing".

### 8. Missing tests to add now
Each: proposed test name in `function_returns_expected_on_condition` form, input class covered, specific bug it would catch, and the test as code or specification. Grouped by function under test. The `reliability` sub-agent's challenge-tests output feeds this section directly.

### 9. What's good
Up to 5 specific, transferable patterns worth keeping, each one sentence with a file reference. No general praise. Skip the section entirely if nothing specific qualifies.

### 10. Commands to re-verify
- Commands the reviewer ran (re-run to confirm they still pass).
- New commands or test invocations the review introduced.

### Author response convention
Address each finding by its identifier (e.g., "B2", "M5") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions from section 4 first.

---

Be direct. If something is wrong, say so plainly and show the fix. Vague praise and vague criticism are equally useless.

## Saving the report

### Directory and filename

Save to the project's `doc/devel/reports/reviews/` directory:

```
doc/devel/reports/reviews/<module-slug>_<YYYY-MM-DD>.md
```

Examples:

- `doc/devel/reports/reviews/gvcf_parser_2026-04-13.md`
- `doc/devel/reports/reviews/genotype_merging_2026-04-13.md`
- `doc/devel/reports/reviews/pr-142_2026-04-13.md`

If a review for the same scope and date already exists, append `_v<N>`:

- `doc/devel/reports/reviews/gvcf_parser_2026-04-13_v2.md`

### Document header

```markdown
# Code Review: <module-slug>
**Date:** <YYYY-MM-DD>
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** <one-line description>
**Status:** <Approve / Approve-with-changes / Request-changes>

---
```

The body is sections 1–10 of *Output format* above, in order, with verbatim headings.

### File links inside findings

References to source files use relative Markdown links from the `doc/devel/reports/reviews/` directory:

- Single line: `[file.rs](../../../../path/file.rs#L123)`
- Range: `[file.rs](../../../../path/file.rs#L123-L456)`

Display text is the path (no backticks).

### Pre-save checklist

- [ ] Every Blocker finding has High confidence; lower-confidence Blocker-class issues are filed at Major with a verification step.
- [ ] Every Major and Minor finding has a concrete fix (code, test, or specific question).
- [ ] Every cited file:line was actually read (no invented locations).
- [ ] Test recommendations use `function_returns_expected_on_condition` naming.
- [ ] "What's good" has 3–5 specific patterns or is omitted.
- [ ] File paths are repo-relative (`src/foo.rs`, not `/home/...`).
- [ ] Severity codes (B1, M1, Mi1) are consistent and dense (no gaps).
- [ ] Open questions are numbered and referenced from the findings they affect.
- [ ] Per-category files in `tmp/review_<date>_<slug>/` are left in place as an audit trail.
- [ ] `PROJECT_STATUS.md` updated (per *Project status protocol*).

## Project status protocol

The project tracks the lifecycle of every feature in `PROJECT_STATUS.md` at the project root. It is a navigation aid, not a source of truth — use it to find the relevant spec, plan, and prior reports for the in-scope feature, then verify against current code as usual.

**At task start.** Read `PROJECT_STATUS.md`. The immutable "About this project" paragraph (delimited by `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` HTML comments) gives the design context and points at the authoritative spec; "Current focus" confirms the project's direction and last-completed work; the per-feature blocks point at the plan, prior impl reports, and reviews for the in-scope feature.

**At task end** (after the review report has been saved): update only the in-scope feature's block in `PROJECT_STATUS.md`.

- Append a link to the new review under `Latest review:` (or `Latest reviews:` if multiple). Prefer to replace the previous link rather than accumulate a long list — `git log` and `ls doc/devel/reports/reviews/` carry chronology.
- Update `Status:` to `reviewed` (or keep `fixes-applied` / `shipped` if the review found nothing material; explain in one trailing word: `shipped (re-reviewed 2026-MM-DD)`).
- Add any new `Open:` items the review surfaced; do not close existing ones (a review surfaces, it does not resolve — that is the apply-fixes skill's job).
- Refresh **Current focus** — rewrite `Last completed task` to name this review and link the report. Touch `Next task` only if the human PM has not already set one; otherwise leave it alone, optionally appending `(suggested follow-up: apply fixes)`.

**Status vocabulary:** `planned` / `in-flight` / `implemented` / `reviewed` / `fixes-applied` / `shipped` / `superseded`. After a code-review run, the typical new status is `reviewed`.

**Hard rules.**

- Do not edit the **About this project** paragraph or anything between the `ABOUT-PARAGRAPH-START` / `ABOUT-PARAGRAPH-END` comments.
- Do not modify another feature's block.
- Do not summarize the review's findings inside the block — the block is a list of pointers; findings live in the saved report.
- If the in-scope feature has no block yet, create one using the format of existing blocks.
- If `PROJECT_STATUS.md` and the current code disagree, trust the code; the status file is stale and should be updated, not relied on.

## Reusable prompt template

Use this to invoke the skill consistently. Fill in the Context block; the rest defers to the skill body.

> Perform a Rust code review per the **rust-code-review** skill. Follow its principles, procedure, severity rubric, and output format in full — do not abbreviate or skip sections.
>
> **Context**
> - **Scope:** <files / module / PR diff / branch comparison>
> - **Domain intent:** <one-paragraph description of what this code is meant to do>
> - **Audience:** <internal service / public library / CLI / embedded>
> - **Constraints:** <performance budgets, MSRV, `no_std`, target platforms, deadline pressure>
> - **Out of scope:** <legacy modules being deleted, generated code, vendored deps>
> - **Prior review history:** <previously reviewed? known tracked issues?>
>
> **Anti-hallucination contract.** Quote tool output verbatim. If a command was not run, list it under "Commands not run". If a file or line cannot be located, say so. Never invent file paths, line numbers, error messages, clippy warnings, test results, or behavior. Findings without verifiable evidence are labeled "Needs verification" per the skill.
>
> **Reminders of the most-violated rules** (not a substitute for the per-category checklists):
> 1. Reliability first — verify behavior, not style.
> 2. Errors must never pass silently.
> 3. Tests cover edge cases and every error path.
> 4. Names are precise, domain-relevant, verb-based for functions.
> 5. Defaults are visible at the call site, in docs, and at runtime.
> 6. Code smells get concrete refactors, not vague complaints.
> 7. Make the compiler flag refactors — no `..Default::default()`, exhaustive destructures.
