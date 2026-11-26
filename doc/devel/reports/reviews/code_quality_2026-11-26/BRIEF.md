# Code-quality map — brief shared by every mapping agent (2026-11-26)

You are one of several agents mapping one area of the `pvc` variant caller (Rust) for a
code-quality program. This phase is **read-only**: you produce a findings report; later
phases fix things, one area at a time, each fix reviewed.

Work in the worktree `/Users/jose/devel/pop_var_caller-code-quality`. Read its `CLAUDE.md`,
`doc/devel/specs/design_principles.md`, and the arch docs under `doc/devel/ng/arch/` that
cover your area (note: `module_layout.md` is stale — it describes a `src/ng/` tree that
no longer exists). Read `ai/skills/clear-technical-writing/SKILL.md` before writing.

## The owner's goals
1. **Architecture**: more isolated modules — *where that does not degrade performance or
   make the code harder to think about*. Wrong-way dependencies (library code importing
   `crate::run` or `crate::cli`, mutual dependencies between modules) are prime targets.
2. **Duplication**: find it; propose what to merge.
3. **Bugs**: find them.
4. **Readability**: code easier to read, understand, reason about.

## The owner's rule — apply it to every finding
**There are no absolute rules; weigh pros and cons.** More modularity can make code
*harder* to reason about; deduplicating two similar functions can couple things that
should evolve apart; splitting a big file can scatter one idea across five. Every
proposal must state its downside, and say "leave it" when the downside wins. A short
list of high-value findings beats a long list of mechanical ones.

The program can do autonomously anything that does not change output, slow the caller,
or raise memory use. Flag clearly any proposal that would do one of those.

## Constraints
- **Do not build, run cargo, or run tests** — other agents share the machine. Reason
  from code. (A separate agent is building an output-identity check.)
- Write only your one report file (path given in your prompt). Change no source.
- Test code (`#[cfg(test)]` blocks, `tests.rs`) is in scope only for duplicated test
  setup and for tests that do not test what they claim; do not review it line by line.
- Inside the crate, public interfaces may change freely; the owner wants them minimal.
  The command line and output files must not change (except by an approved bug fix).

## For bugs specifically
A bug finding needs a **concrete failure scenario**: inputs/state → wrong output, panic,
or silent data loss. Rank by consequence to calls (wrong genotype/QUAL/missing record >
crash on rare input > wrong log message). Mark each **suspected** (reasoned) or
**likely** (you traced the path end to end). Think about the edges the caller must
survive: one sample, thousands of samples, 3x and 300x depth, empty regions, contig
ends, zero reads, multi-allelic sites, overlapping records, integer overflow, float
NaN/inf, thread-count dependence, ordering dependence.

## Report format (Markdown)
Start with a 5–10 line summary: the top findings in priority order, in plain words.
Then one section per category (Bugs, Architecture, Duplication, Readability). Each
finding:

- **Title** — one line, plain English.
- **Where** — `path:line` (verify every citation by reading it; wrong lines waste the fixer's time).
- **What** — the problem, with its size (how many lines/call sites/copies).
- **Proposal** — what to do.
- **Downside / why it might be wrong** — honest.
- **Effect** — output change? speed/memory risk? (none / possible / yes)
- **Cost** — S (<1h), M (half day), L (more).
- **Confidence** — high / medium / low.

End with "Looked at, left alone" — things that look odd but are fine, one line each, so
the next agent does not re-litigate them.

Your final chat reply to the manager: the report path and the summary, nothing more.
