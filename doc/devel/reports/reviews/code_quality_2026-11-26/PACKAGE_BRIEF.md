# Brief for an agent executing one work package (code-quality program, 2026-11-26)

You execute **one package** from `PLAN.md` (same directory) — its letter-number is in your
prompt. Read, in order: the repo `CLAUDE.md`; `PLAN.md` "How to use this plan" and your
package's row; the report sections your row cites (they hold the evidence and the exact
locations); `ai/skills/clear-technical-writing/SKILL.md` before writing commit messages or
doc comments.

## Where to work

Create your own worktree from the integration branch and work only there:

```
git -C /Users/jose/devel/pop_var_caller-code-quality worktree add -b cq-<pkg> /Users/jose/devel/pop_var_caller-cq-<pkg> code-quality
```

(`<pkg>` lowercase, e.g. `cq-f1`.) Use absolute paths in every command — the shell's working
directory resets between calls. Scratch goes in your worktree's `tmp/`, never `/tmp` or a
harness scratchpad. Python via `uv run`, never `python3`. Do not push. Do not touch other
worktrees. Do not merge into `code-quality` — the manager does that.

## The owner's rules

- **No absolute rules: weigh pros and cons.** If, once in the code, an item turns out not worth
  doing (the downside wins, or the finding was wrong), skip it and say why. That is a good outcome.
- You act alone only on changes that **alter no output, slow nothing, and raise no memory**.
  If an item turns out to change output on valid input, or the gate flags time/memory, stop that
  item, revert it, and report it with the evidence — it becomes an owner decision.
- Stay inside your package's files. If a fix needs a file another package owns, report it
  instead of editing.

## How to build and test (shared machine: other agents build at the same time)

- Native cargo, your own target dir: `CARGO_TARGET_DIR=<your worktree>/target-host`,
  `CARGO_BUILD_JOBS=4`.
- Light dev loop: run the tests of the modules you touched (`cargo test --release --lib <module path>`
  and any integration test that covers them), plus `cargo clippy --all-targets` and
  `cargo fmt --check` *on the files you touched only* (rustfmt those files, never
  `cargo fmt` on the whole crate).
- **Bug fix commits**: one bug per commit; the commit adds a test that fails before the fix
  (show that it failed — run it on the unfixed code first) and passes after.
- **Every commit's output is checked** with the gate in identity mode:
  `<your worktree>/scripts/code_quality_oracle.sh --identity` (~2 min, outputs only).
  It must exit 0. Do NOT run the timing mode (no flag) — timing is measured serially by the
  manager after merging, because parallel builds make timings meaningless. If your package
  touches a hot loop, say so in your reply so the manager times it specifically.
- Things the gate cannot see (log text, errors, `--help`, CRAM, long reads, many read groups,
  thousands of samples) need their own test or a before/after diff in the commit message.

## Commits

Small, one purpose each; bug fixes separate from refactors. Message: what changed and why, in
plain words; end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Final reply to the manager (short)

Branch name and commit list (hash + subject); for each plan item: done / skipped (why) /
stopped for the owner (evidence); gate result; anything the manager must time or decide;
any finding that belongs to another package.
