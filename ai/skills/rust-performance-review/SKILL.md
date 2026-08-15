---
name: rust-performance-review
description: Use this skill whenever the user asks for a performance review, hot-path audit, profiling guidance, or to "find places to make this faster" in Rust code. Trigger on phrases like "look for performance improvements", "is this on the hot path", "where is this allocating", "review for cache locality", "is there lock contention here", or when the user shares a flamegraph / criterion result and asks for next steps.
---

# Rust Performance Review

You are scanning Rust code for performance-improvement candidates. You are **not** rewriting the code, and you are **not** asserting wins; you are surfacing places where a measurable gain is plausible and proposing the experiment that would confirm or refute it.

This work has a hostile prior: most "obvious" optimizations either do not matter (cold code), do not help (the compiler already handles it), or trade real complexity for an imagined win. Be skeptical. Each candidate must answer: where on the call graph does this run? what would we measure? how much complexity does the fix add?

The review is split across focused per-category checklists in `ai/skills/rust-performance-review/performance_review/`. **You are the orchestrator**: you triage which categories apply to the scope, dispatch one sub-agent per category in parallel — **each in its own git worktree** (step 5) — then synthesize their findings into a single report. Per-category rules are not duplicated in this file — read each category file when dispatching the corresponding sub-agent.

## Review principles (must always hold)

- **Profile first; pattern-match second.** A finding without a profile, benchmark, or strong call-graph evidence that the code is on a hot path is filed at low priority and paired with the measurement that would promote it. Do not propose a rewrite of code you cannot show is hot.
- **Every candidate has a measurement plan.** State the experiment: which benchmark, which profiler, which metric (wall time / allocations / cache misses / lock-wait time / syscalls), and what threshold makes the fix worth merging. "It will be faster" is not a finding.
- **Complexity is a cost.** Every fix names the complexity it introduces (extra type, lifetime gymnastics, `unsafe`, build-config knob, dependency) and weighs it against the expected gain. A fix that doubles maintenance for a 2% wall-time win is a bad trade in critical lab code.
- **Hot-path discipline.** Optimizations in code that runs once at startup, in a CLI flag handler, or in error-handling paths are noise. Be explicit about the call-frequency assumption and downgrade severity when call frequency is unverified.
- **Correctness first.** Never recommend an optimization that weakens invariants, introduces `unsafe`, or relaxes atomic ordering without a separate, justified safety review. Performance findings that touch correctness boundaries are flagged as such and held to the evidence bar of a correctness review.
- **One change per measurement.** Bundling allocator switch + LTO + a code refactor in the same PR produces an unreadable result. Each PR names the single hypothesis being tested.

The severity rubric and per-finding format are defined in `ai/skills/rust-performance-review/performance_review/_finding_format.md`. Read it once at the start of every review — both you (for synthesis) and every sub-agent you dispatch will follow it.

## Review procedure

You run each numbered step once per review. Only step 5 fans out into parallel sub-agents.

### 0. Orient on prior work for this scope

Before reviewing, find what already exists for the code in scope, so the review starts where the last one left off instead of re-deriving it. Look for: a prior performance review of the same module (`doc/devel/reports/reviews/perf_*.md` — the filename prefix `perf_` separates them from correctness reviews), the spec and architecture docs (`doc/devel/specs/`), the implementation plan (`doc/devel/implementation_plans/`), and any implementation report (`doc/devel/reports/implementations/`). A prior perf review's open findings and its measurement plan are the most useful of these: they tell you which experiments were already run and what they showed.

### 1. Establish scope and call frequency

Determine whether the review covers a full crate, a module, a PR diff, or a single function under benchmark. State the scope. For each in-scope file, identify which functions are believed to be on the hot path and on what evidence: a profile, a benchmark, a call-graph argument, or "no evidence yet". Code with no hot-path evidence is reviewed at lower priority — flag this as a constraint passed into every sub-agent prompt.

### 2. Inventory existing measurement

List every existing benchmark (`benches/`), profile artifact (`flamegraph.svg`, `*.profraw`, `dhat-*.json`, `samply.json.gz`), and any quoted measurement the user has provided. Note what is missing — which functions have no benchmark, no profile, or only synthetic timings. **If no measurement exists for the in-scope code, the first deliverable is a measurement plan, not code rewrites** — say so explicitly in the verdict.

If real benchmark or profile output is available, quote it verbatim. The verbatim output is passed into each sub-agent's prompt at step 5 so they do not re-run.

**If a sampling profile is impossible to collect in the user's environment, raise the alarm before moving on.** This is not a "we'll work around it with other tools" situation — sampling profile access is a load-bearing input to a useful review, and proceeding without it significantly degrades the review's actionable output.

**Read `ai/skills/rust-performance-review/performance_review/profiling_environment.md` in this step, before writing any measurement plan.** It holds the machine-specific inventory this step depends on: which profiling tools are verified to work on the macOS host and in the Linux dev container (and their working invocations), what this machine cannot measure (the PMU is not virtualized, so hardware counters — cache-miss, branch-miss, `perf c2c`, `perf sched` — are unavailable on host and container alike, with the software stand-ins listed there), the per-OS diagnosis for a blocked sampling profiler, fixture locations, and the toolchain pin. Every sub-agent dispatched at step 5 is pointed at the same file — a measurement plan that names a tool the machine cannot run is a plan the review cannot execute.

**Other profilers are not substitutes for a sampling profile.** Each of them answers one narrow question:

- DHAT — allocations only. A function can be 30 % of CPU and 0 % of allocations.
- `cargo asm` — codegen inspection. Tells you whether bounds checks elided or autovec happened, but not whether the affected line matters.
- valgrind / callgrind — deterministic instruction counts. Useful for ranking changes against each other but ~10–50× slower than native; the instruction count is *not* the same as the CPU-cycle cost (cache misses, branch mispredicts, FMA throughput all invisible).
- Criterion wall-clock — cross-run variance routinely exceeds the effect size of individual changes: system load can shift the median 20–25 % between runs even when the within-run confidence interval reads 1 % (see `methodology.md` on cross-commit comparisons and the revert experiment).

None of these tells you which line of code holds 30 % of the CPU. Only a sampling profile does. A review without one will produce a long list of pattern-matched candidates of which most will show **no measurable gain when applied**, because the underlying site was not in the top 20 % of self-time. This failure mode is documented end-to-end in `doc/devel/reports/reviews/perf_baq_2026-05-12.md` (round one applied seven Likely findings against pattern-match alone; the criterion bench could not detect a net change; round two used a real profile and got measurable wins on the same code base).

**In a pipeline, pair the sampling profile with one off-CPU signal before trusting its ranking.** A sampling profile ranks on-CPU self-time, but in a producer→workers→writer pipeline the wall clock is set by whichever stage the others wait on — and a stage can wait off-CPU (blocked channel send/recv, lock wait, I/O), where the sampler records nothing. When the scope is pipeline-shaped, collect one off-CPU measurement alongside the CPU profile: channel queue depth and blocked-send vs blocked-recv time (`hotpath` instrumentation or a few manual timers), or a causal profile (`coz`) — invocations in `profiling_environment.md`. A queue that sits full in front of one stage names the bottleneck more directly than any self-time ranking; a flamegraph of a blocked pipeline shows where CPU goes, not where wall time goes.

When sampling is blocked, **tell the user explicitly**, in this order:

1. Name what is blocked and why (kernel paranoid setting, missing capabilities, sandboxed runner — diagnose precisely).
2. Name the concrete fix to apply (the `sysctl` line, the container flag, "run on bare metal", etc.).
3. State the cost of *not* fixing it: the review's actionable output is significantly degraded; the candidate list will not be rankable by hot-path evidence and many candidates that look promising will fall through.
4. Ask whether they want to fix the access issue first, or proceed at lower-quality output.

Do not silently route around this with "we'll use DHAT instead" or "we'll use cargo asm and infer." If the user accepts proceeding without sampling, downgrade every code-level finding's severity to **Likely** (no **Hot-path** is reachable without profile evidence per the rubric) and call out the constraint at the top of the report.

### 3. Determine intent and targets

Before judging whether the code is fast enough, establish what it is meant to do — its purpose, its expected input sizes, its latency or throughput target, its target hardware. Performance review is meaningless without these numbers. If the targets are not stated, ask before continuing or file the gap as a Note finding.

Write a one-paragraph "performance intent" summary; it is passed into each sub-agent's prompt.

### 4. Triage categories

Decide which per-category checklists apply. Each lives at `ai/skills/rust-performance-review/performance_review/<category>.md`.

| Category | Apply when |
|---|---|
| `methodology` | Always. Establishes that profiling, benchmarking, and `Cargo.toml` are sound before any code-level finding is acted on. |
| `allocations` | Code on the hot path constructs `Vec` / `String` / `Box` / `Arc` / map types, clones owned data, calls `format!` / `to_owned` / `to_string`, or has unbounded buffer growth. |
| `data_layout` | Hot path iterates over collections of structs, holds shared atomic state across threads, defines a `pub struct` with several fields, or there is suspicion of cache-miss pressure or false sharing. |
| `concurrency` | Code uses `Arc`, `Mutex`, `RwLock`, atomics, channels, `rayon`, `tokio`, `async fn`, or `spawn`. Skip otherwise. |
| `hot_loops` | Hot path contains tight numeric or byte-processing loops — including float reductions (sums, dot products, likelihood accumulations) and filters branching on data-dependent predicates — slice indexing, iterator chains, generic dispatch, `format!`, or anywhere autovectorization or branch layout could plausibly matter. |
| `io_and_syscalls` | Code performs file or socket I/O, reads/writes large data, or makes per-record syscalls. |

When in doubt, dispatch — a sub-agent that finds nothing applicable writes `No findings.` and is cheap.

### 5. Dispatch sub-agents in parallel — **each in its own worktree**

Create the scratch directory: `tmp/perf_review_<YYYY-MM-DD>_<scope-slug>/` (append `_v<N>` if it already exists). The slug is a short kebab-case identifier of the reviewed module or PR (e.g. `gvcf_parser`, `pr-142`).

Project rule: scratch space is project-local `tmp/`, never `/tmp`. Add `tmp/` to `.gitignore` if it is not already covered by the existing target ignores.

For each selected category, dispatch a `general-purpose` sub-agent **in parallel** — issue a single message with multiple Agent tool calls — and pass **`isolation: "worktree"`** on every one.

#### Why the isolation is mandatory here in particular

This skill's whole currency is measurement, and **a measurement taken on a shared checkout is not a measurement**. Two agents benchmarking at once contend for the same target directory and the same machine; one rebuilding while another times is an unattributable number, and an ablation is worse — it needs the tree to hold *one* change while a baseline and a variant are timed in turn.

A code review of a shared tree has already been observed losing to this: nine agents, edits overwritten mid-experiment, `src/` reverted under an agent three times, one agent's baseline failing on another's marker, five of nine retreating into private worktrees late and after wasted work, and every mutation result needing serial re-verification before it could be believed. A perf review is more exposed, not less, because its results are continuous — interference does not show up as a failed build, it shows up as a plausible wrong number.

Two further rules follow, and they are this skill's not the code review's:

- **Ablations and baselines belong to one agent inside one worktree.** Never split "measure the baseline" and "measure the variant" across agents; they will not be comparable.
- **Say whether the numbers are comparable *across* agents.** Separate worktrees mean separate target directories and concurrent load. Timings from different agents are the same order of magnitude, not the same experiment — at synthesis, treat a cross-agent comparison as needing a re-run by one agent before it enters the report.

#### The worktree starts on `main`, not on your branch — every prompt needs a step 0

**`isolation: "worktree"` checks the agent's tree out at `main`.** Probed during a review on a
feature branch: the agent's `HEAD` came back exactly `git rev-parse main`, and the module under
review was missing two of its files there. **The code under review is absent from every agent's
tree**, so an agent handed its scope by absolute path silently works in *your* checkout instead —
the one place those files exist. On a perf review that is worse than on a code review: several
agents benchmarking in one tree produce plausible wrong numbers rather than a failed build.

So **every agent prompt begins with a step 0**: `git checkout --detach <sha>` in its own worktree
root, then `ls` two files that exist only on the branch, and **stop and report** if either is
missing. Verified: the branch-only files were absent before and present after, the tree stayed
clean, and no fetch was needed — the worktrees share one object store. An **unchanged** worktree is
also auto-cleaned, so a resumed agent can land in a real checkout, which is why the check is a hard
stop rather than a warning.

#### What isolation changes in the prompt

1. **Where it builds and benchmarks — its own copy of the build script, by that copy's absolute
   path.** The agent is *in* its worktree, and `scripts/dev.sh` derives `PROJECT_DIR` from the
   script's own location: **the agent's copy mounts and builds the agent's tree**, while
   `<main-checkout>/scripts/dev.sh` builds the main checkout and ignores the caller's directory —
   which on this skill means benchmarking code that is not the code under review. Instruct the
   agent by its own worktree's absolute path. (An earlier note claiming `dev.sh` "cannot build a
   worktree" described an invocation of the *main* copy; it is wrong, and it misdirected a
   fan-out.)
2. **Where its findings go.** Its worktree is temporary and is cleaned up, so the findings file goes to an **absolute path in the main checkout's** `tmp/perf_review_<date>_<slug>/`. Each agent writes its own file.
3. **What it may now do.** It has its own tree: it may apply a candidate fix and measure it rather than only proposing a measurement plan. Say so — a measured fix outranks a plan, and this is what buys it.

Do **not** ask isolated agents to leave the tree clean or revert their experiments — that is a shared-checkout instruction, and here it wastes effort on a tree about to be discarded.

Each sub-agent prompt:

> Run the **<category>** checklist on the following Rust performance-review scope.
>
> **Scope:** <full crate / PR diff / module / single function>
> **Performance intent and targets:** <one paragraph from step 3, including throughput/latency targets and target hardware>
> **In-scope files (full paths):** <list>
> **Hot-path evidence:** <verbatim profile output / benchmark numbers / "none — pattern-match only">
> **Out of scope:** <list, with reasons>
>
> **Instructions:**
> 0. **Re-point your worktree at the commit under review.** It was created from `main` and does
>    not contain this branch's code. Run `git checkout --detach <sha>` in your own worktree root,
>    then confirm <two branch-only files> exist. If either is missing, **stop and report** rather
>    than measuring another checkout.
> 1. Read `ai/skills/rust-performance-review/performance_review/<category>.md` for the rules to apply.
> 2. Read `ai/skills/rust-performance-review/performance_review/_finding_format.md` for the severity rubric and finding format.
> 3. Read `ai/skills/rust-performance-review/performance_review/profiling_environment.md` for what this machine can and cannot measure. Every measurement plan you write must name only tools that file confirms available — no hardware counters exist here (no `perf stat -e branch-misses`, no `perf c2c`), so use the software stand-ins it lists.
> 4. Read each in-scope file.
> 5. Apply each rule and produce findings in the specified format. For every candidate, propose the **measurement plan** (benchmark or profile that confirms the gain) and the **complexity cost** of the fix. Findings without a measurement plan are downgraded.
> 6. **You are in your own isolated git worktree.** Build, benchmark and mutate freely *there* — nothing you do can collide with another agent. Run the project's tooling **from your own worktree's copy, by its absolute path** (`<your-worktree>/scripts/dev.sh …`): the script builds whichever tree it lives in, so the main checkout's copy would measure the wrong code. Where you can, go past proposing a measurement plan and **run it**: apply the candidate, measure baseline against variant *in this one worktree*, and report both numbers. A measured fix outranks a proposed one. Keep every ablation inside your own tree; a baseline timed in one worktree and a variant in another are not comparable.
> 7. Write findings to the **absolute path** `<main-checkout>/tmp/perf_review_<date>_<slug>/<category>.md` — outside your worktree, which is temporary and will be discarded. If no findings apply, write only the line `No findings.`
> 8. Do not invent profile output, benchmark numbers, or call frequencies. Cite only what was provided in the prompt, what you read in the source, or what you measured yourself — and say which.
> 9. Stay within the category. Issues that belong elsewhere go under a `## Cross-category observations` heading at the bottom of your file.

Substitute `<category>` and the scope fields for each dispatch. Do **not** assign severity codes (H1, L1, …) inside sub-agents — that happens at synthesis.

### 6. Collect findings

When all sub-agents complete, read each `tmp/perf_review_<date>_<slug>/<category>.md`. Tally findings, note any cross-category observations, and decide whether each cross-category note becomes its own finding or merely informs synthesis.

Promote findings whose measurement plans converge (multiple categories agree the same site is hot). Demote findings that turn out to be cold-path on closer inspection. If a sub-agent appears to have skipped its scope or produced unverifiable findings, redispatch that one category with an explicit instruction to fix the gap.

### 7. Synthesize the unified report

Compose the report using the *Output format* below. Verdict, measurement plan, and "What's already good" need the full picture and are produced by you. Assign severity codes during synthesis: `H1, H2, …` for Hot-path; `L1, L2, …` for Likely; `S1, S2, …` for Speculative; Notes stay grouped without numbering. Save to `doc/devel/reports/reviews/perf_<module-slug>_<YYYY-MM-DD>.md` per the saving conventions below. Leave the per-category files in `tmp/` as an audit trail.

## Output format

Produce the synthesized report in the following order. Use the section headings verbatim so the format is machine-readable.

### 1. Scope and constraints
- What was reviewed: full crate / module / PR diff / single function.
- Reviewed against: commit hash, branch, or "as-provided".
- Throughput / latency targets, expected input sizes, target hardware.
- Hot-path evidence available (profile, benchmark, or "none — pattern-match only").
- In-scope files (list).
- Deliberately out of scope (list, with reason).
- Categories dispatched (list, each with a one-line reason).

### 2. Verdict
One of:
- **Profile first** — there is not enough hot-path evidence to recommend code changes. Section 3 lists the measurements to take.
- **Run experiments** — candidates exist; their priority and order are listed in section 5.
- **Apply the listed wins** — at least one candidate is well-evidenced (matching profile, plausible mechanism, contained complexity); apply with the proposed measurements as gates.

### 3. Measurement plan
The benchmarks and profiles to add or run, in the order they unblock other findings. Each entry: command, expected output, what threshold answers what question. If the verdict is *Profile first*, this is the primary deliverable.

### 4. Build / toolchain configuration
LTO, codegen-units, panic, opt-level, debug, allocator, target-cpu — anything in `Cargo.toml` or `.cargo/config.toml` that should change before code-level work. Driven by the `methodology` sub-agent's findings.

### 5. Code-level findings
Grouped by severity (**Hot-path** → **Likely** → **Speculative** → **Note**). Within each severity, ordered by confidence (High → Low), then by file. Each finding follows the format defined in `ai/skills/rust-performance-review/performance_review/_finding_format.md`, with the severity code (H1, L1, …) prepended to the title — e.g. `H1: src/parser.rs:42 — Title`.

### 6. Out-of-scope observations
Performance smells in untouched code, surfaced but not blocking. Each: file, brief description, suggested follow-up (separate PR or issue).

### 7. What's already good
Up to 3 specific, transferable patterns the code is already getting right (e.g., "uses `with_capacity` everywhere it knows the size", "shards by sample to avoid `Mutex` contention"). One sentence each, with file references. No general praise. Skip the section entirely if nothing specific qualifies.

### Author response convention
Address each finding by its identifier (e.g., "H1", "L2") with one of: `applied in <commit>` / `experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. The "experiment shows no gain" path is expected and welcome — that is what the measurement plan is for.

---

Be direct. If something is plausibly hot, say so plainly and propose the experiment. Vague "could be faster" and vague "looks fine" are equally useless.

## Saving the report

### Directory and filename

Save to the project's `doc/devel/reports/reviews/` directory:

```
doc/devel/reports/reviews/perf_<module-slug>_<YYYY-MM-DD>.md
```

Examples:

- `doc/devel/reports/reviews/perf_gvcf_parser_2026-05-10.md`
- `doc/devel/reports/reviews/perf_pipeline_2026-05-10.md`
- `doc/devel/reports/reviews/perf_pr-142_2026-05-10.md`

If a review for the same scope and date already exists, append `_v<N>`:

- `doc/devel/reports/reviews/perf_gvcf_parser_2026-05-10_v2.md`

### Document header

```markdown
# Performance Review: <module-slug>
**Date:** <YYYY-MM-DD>
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** <one-line description>
**Verdict:** <Profile first / Run experiments / Apply the listed wins>
**Hot-path evidence:** <profile / benchmark / pattern-match only>

---
```

The body is sections 1–7 of *Output format* above, in order, with verbatim headings.

### File links inside findings

References to source files use relative Markdown links from the `doc/devel/reports/reviews/` directory:

- Single line: `[file.rs](../../../../path/file.rs#L123)`
- Range: `[file.rs](../../../../path/file.rs#L123-L456)`

Display text is the path (no backticks).

### Pre-save checklist

- [ ] Every Hot-path finding has High confidence and quotes the profile / benchmark output that names the site.
- [ ] Every Likely finding has a measurement plan that would confirm or refute the gain.
- [ ] Every finding has a complexity cost named honestly.
- [ ] Every cited file:line was actually read (no invented locations).
- [ ] No fabricated percentages or speedup multipliers anywhere in the report.
- [ ] Cold-path findings are downgraded to Note unless the orchestrator marked them in scope.
- [ ] Severity codes (H1, L1, S1) are consistent and dense (no gaps).
- [ ] Per-category files in `tmp/perf_review_<date>_<slug>/` are left in place as an audit trail.
- [ ] Build configuration findings (section 4) are separated from code-level findings (section 5).

## Reusable prompt template

Use this to invoke the skill consistently. Fill in the Context block; the rest defers to the skill body.

> Perform a Rust performance review per the **rust-performance-review** skill. Follow its principles, procedure, severity rubric, and output format in full — do not abbreviate or skip sections.
>
> **Context**
> - **Scope:** <files / module / PR diff / branch>
> - **Performance intent and targets:** <what the code does, on what input sizes, with what latency / throughput target, on what hardware>
> - **Hot-path evidence:** <flamegraph path / criterion bench output / DHAT report / "none — pattern-match only">
> - **Constraints:** <MSRV, `no_std`, target platforms, deadline pressure>
> - **Out of scope:** <legacy modules being deleted, vendored deps>
> - **Prior review history:** <previously reviewed? known tracked issues?>
>
> **Anti-hallucination contract.** Quote tool output verbatim. If a benchmark or profile was not run / read, say so under "Hot-path evidence" and downgrade findings accordingly. Never invent measurements, percentages, file paths, or line numbers.
>
> **Reminders of the most-violated rules** (not a substitute for the per-category checklists):
> 1. Profile first; pattern-match second.
> 2. Every finding has a measurement plan and a complexity cost.
> 3. Cold code is not the hot path; downgrade findings filed against it.
> 4. Correctness wins over speed. Route correctness-adjacent findings through a separate review.
> 5. One change per measurement — do not bundle allocator + LTO + refactor.
> 6. **If the environment blocks sampling profilers (`perf record` / `cargo flamegraph` / `samply` / Instruments / `dtrace`), halt and tell the user before continuing.** This is not a "we'll use DHAT instead" situation — DHAT, `cargo asm`, valgrind, and wall-clock criterion each answer one narrow question and **none** of them tells you which line owns the CPU. Without a sampling profile, most pattern-matched Likely findings will show no measurable gain when applied. The fix depends on OS (Linux: `sudo sysctl kernel.perf_event_paranoid=2`; macOS: usually no fix needed — `samply` / Instruments work without elevated perms); see step 2 and `performance_review/profiling_environment.md` for the per-OS diagnosis and this machine's verified inventory.

---

## Sources

The two articles that motivated this skill:

- *Inside Rust's std and parking_lot Mutexes: Who Wins?* — https://blog.cuongle.dev/p/inside-rusts-std-and-parking-lot-mutexes-who-win
- *About memory pressure, lock contention and data-oriented design* — https://mnt.io/articles/about-memory-pressure-lock-contention-and-data-oriented-design/

Background reading the per-category rules draw on:

- *The Rust Performance Book* (Nicholas Nethercote) — https://nnethercote.github.io/perf-book/
- *How to avoid bounds checks in Rust (without unsafe!)* — https://shnatsel.medium.com/how-to-avoid-bounds-checks-in-rust-without-unsafe-f65e618b4c1e
- *Why my Rust benchmarks were wrong* (Guillaume Endignoux on `black_box`) — https://gendignoux.com/blog/2022/01/31/rust-benchmarks.html
- *Optimization adventures: making a parallel Rust workload faster with (or without) Rayon* — https://gendignoux.com/blog/2024/11/18/rust-rayon-optimized.html
- `crossbeam_utils::CachePadded` docs — https://docs.rs/crossbeam-utils/latest/crossbeam_utils/struct.CachePadded.html
- `bumpalo` — https://github.com/fitzgen/bumpalo
- `dashmap` — https://github.com/xacrimon/dashmap
- `rustc-hash` (FxHash) — https://github.com/rust-lang/rustc-hash
- `mimalloc` — https://github.com/microsoft/mimalloc
- `flamegraph-rs` — https://github.com/flamegraph-rs/flamegraph
- `criterion.rs` — https://github.com/bheisler/criterion.rs
- `hotpath` (off-CPU channel/lock/async instrumentation; no PMU required) — https://hotpath.rs, *Profiling async Rust* — https://hotpath.rs/blog/profiling-async-rust, and *A guide to profiling Rust* (deterministic signals over timing; off-CPU layers before CPU sampling in pipeline-shaped services) — https://hotpath.rs/blog/profiling-rust-guide
- *Speeding up float calculations in Rust* (algebraic float ops; pairwise summation) — https://pythonspeed.com/articles/faster-float-math-rust/
- *Branchless programming in Rust* (misprediction cost; the always-write/conditionally-advance filter) — https://www.greyblake.com/blog/branchless-rust/
- `fearless_simd` (safe SIMD with runtime dispatch) — https://linebender.org/blog/fearless-simd-0-7/
