# Project operating guidelines

These rules apply to every Claude Code session in this repo. They exist
so the assistant can work autonomously without permission prompts on
trivial operations and without accidentally writing outside the
project.

## Container vs. host: where commands run

The project ships a development container (`./scripts/dev.sh`) that
mounts the project tree at the same absolute path inside the
container. **All write-side build tooling runs in the container**, so
nothing you do can touch the host outside this directory.

`$HOME/genomes` is also mounted, **read-only**, whenever it exists —
benchmark configs point reference paths there (e.g. tomato1's
`S_lycopersicum_chromosomes.4.00.fa`), so runs that need a reference
work without extra setup. Read-only keeps the guarantee that matters:
the container can only *write* inside the project. Override the
location with `DEV_GENOMES_DIR`; mount one further path with
`DEV_EXTRA_MOUNT`.

The wrapper picks a runtime automatically: **podman** if present
(Linux dev box), otherwise Apple's **`container`** CLI on macOS
(https://github.com/apple/container).

**Some machines have neither, and `rick` is one of them** — the Linux box the
archive-scale walks run on has no container runtime installed. There
`./scripts/dev.sh` refuses to run rather than falling through to the host,
because a silent fallback would drop the only guarantee it offers. On such a
machine **run `cargo` directly**, and expect two differences:

- the binary lands in `target/release/...`, not `target-container/release/...`
  — the container build points `CARGO_TARGET_DIR` at the second, so the two
  trees can hold different builds of the same example at once. **A script that
  looks for a built binary must check both and take the newer**, which is what
  `scripts/ng_str_library_survey.sh` does;
- `cargo` writes to `~/.cargo` and `~/.rustup` as usual, which the
  containerised run does not.

So **"all write-side build tooling runs in the container" is the rule wherever a
runtime exists, not a fact about every machine.** Do not tell someone working on
`rick` to build with `./scripts/dev.sh`; it will only tell them what is missing.

**One-time macOS setup** (Apple container only):

```
container system start                                      # per-boot
defaults write com.apple.container.defaults build.rosetta -bool false
```

The Rosetta default lets `container build` succeed without installing
Rosetta — fine because we only build arm64 images. The wrapper also
hands Apple container `--memory 16g --cpus 8` because its default VM
size OOM-kills rustc/ld on this project; override with `DEV_MEM=` /
`DEV_CPUS=` env vars if needed. Podman uses host resources directly
and ignores those knobs.

- **Run inside the container** (via `./scripts/dev.sh ...`):
  - Anything that compiles, links, or fetches Rust crates: `cargo
    build`, `cargo test`, `cargo run`, `cargo fmt`, `cargo clippy`,
    `cargo doc`, `cargo fetch`, `cargo search`, `cargo metadata`, etc.
  - Any tool the project depends on that is installed in the
    Containerfile but not on the host.
  - Anything that may write outside the project tree.
- **Run on the host directly** (no container, no permission prompt
  needed because they're allow-listed in `.claude/settings.json`).
  These are granted **only for paths inside this project directory**.
  Do not point them at anywhere else on the host — e.g. don't
  `grep -r ~`, don't `find /etc`, don't
  `cat /Users/jose/somewhere/else`.
  - Read-only inspection of files: `grep`, `rg`, `cat`, `head`,
    `tail`, `less`, `wc`, `ls`, `find` (without `-delete`/`-exec rm`),
    `diff`, `file`, `stat`, `realpath`.
  - Read-only `git` commands: `git status`, `git log`, `git diff`,
    `git show`, `git blame`, `git branch`, `git ls-files`.
  - The dedicated tools (`Read`, `Edit`, `Write`) for file content —
    these are always preferred over shelling out to `cat`/`sed`/etc.

  **Exception — dependency source:** when looking up a dependency's
  API (e.g. noodles), reading from `~/.cargo/registry/src/...` is OK
  because that's just the cargo cache of crates this project depends
  on. It is the only path outside the project tree these read-only
  commands may target. If you find yourself wanting to read anywhere
  else on the host, stop and ask first.

The host allowlist intentionally does *not* include destructive
commands (`rm`, `mv`, write-side `git`) — those still need explicit
approval or should run inside the container.

## Reading dependency source

Cargo extracts crate sources to `~/.cargo/registry/src/index.crates.io-*/<crate>-<version>/`
on the host (and to `/usr/local/cargo/registry/src/...` inside the
container after a build).

**Prefer the Explore agent for API lookups.** When the task is "what
does noodles' Reader API look like?" or "how do I build a CRAM
header with this library?" — delegate to an `Explore` subagent and
ask it to summarise the relevant types, methods, and signatures.
Reasons:

- The agent reads many files across the dependency in parallel and
  returns a compact summary, instead of dumping hundreds of lines of
  third-party source into the main conversation.
- The main context stays focused on the project's own code and
  decisions, not the noodles internals you only needed to look at
  once.
- Multiple lookups (e.g. cram + sam + fasta together) become a
  single agent run rather than a long sequence of `grep`s and reads.

Direct host-side `grep` / `Read` on the registry is still fine for a
one-off check ("does this exact method exist?"). The rule of thumb:
if you'd otherwise read more than two or three dependency files,
spawn an Explore agent instead.

## Scratch space

Do **not** use the system `/tmp`. If you need scratch space for
fixtures, generated outputs, or temporary files, create a project-local
`tmp/` directory (already covered by the surrounding `/target` and
`/target-container` ignores; add `/tmp/` to `.gitignore` if you start
using it). This keeps everything inside the project mount and avoids
leaving state on the host.

**This includes any scratch directory the assistant's own harness offers.**
Claude Code hands each session a scratchpad under the host's
`/private/tmp/...` and tells it to prefer that over `/tmp`; that instruction
is overridden here. Everything — probe scripts, extracted patches, agent
reports, generated fixtures, draft commit messages — goes under this
repository's `tmp/`. Two reasons: a path outside the project mount is
invisible inside the dev container, so anything written there cannot be
handed to `cargo`; and scratch that outlives the session is state left on the
host, which the container sandbox exists to prevent.

## Profiling (samply, flamegraph, perf)

Sampling profilers inside the container use `perf_event_open(2)` on
the host kernel and are gated by `kernel.perf_event_paranoid`. The
container can't relax that; only the host can.

**On the Linux dev box (rootless podman, uid 1000 → uid 1000):**
`--cap-add=SYS_ADMIN` / `--cap-add=PERFMON` on the `podman run`
command line does *not* grant perf access — the syscall reaches the
kernel as the unprivileged invoking user regardless. The only fix is
to lower the sysctl on the host:

```
# user-space sampling (samply on your own binary): paranoid <= 2
# CPU/kernel events (flamegraph, perf record on kernel): paranoid <= 1
sudo sysctl kernel.perf_event_paranoid=1
```

Persistent: drop it in `/etc/sysctl.d/99-perf.conf`. If samply or
flamegraph fails with a permission/`perf_event_open` error inside the
container, this is almost always the cause — don't try to work around
it container-side.

**On macOS (Apple `container`):** containers run inside a lightweight
Linux VM, so `perf_event_open` targets the VM kernel rather than the
host. Linux-only profilers (perf/flamegraph/samply) are best used on
the Linux dev box; for macOS-side profiling use native tools (Xcode
Instruments, `samply` running directly on the host binary).

## Working autonomously

The combination of the container sandbox plus the host read-only
allowlist means most routine work — reading code, searching the
codebase, building and testing inside the container — needs no user
prompts. Use that to keep momentum:

- Don't ask before running `grep` / `cat` / `find` on the host.
- Don't ask before running `cargo *` inside the container.
- Do still pause and ask before destructive operations (force pushes,
  branch deletions, history rewrites, anything that rewrites tracked
  files outside the immediate task).

## Writing for the reader — including in chat

Two skills, and **both are mandatory, not optional reading**:

- `ai/skills/clear-technical-writing/SKILL.md` — all reader-facing prose:
  specs, architecture, doc comments, commit messages, reports.
- `ai/skills/reporting-in-chat/SKILL.md` — **read this before sending any
  reply that reports a result, asks for a decision, or summarises work.**
  It is the conversational counterpart, and it carries a **failure log**:
  append a real before/after every time a reply has to be corrected.

**A reply in chat that explains a result, argues for a design, or
reports an analysis is reader-facing prose.** That is where this rule
keeps being missed: docs and commit messages get the care, and then the
same finding is delivered in chat as compressed shorthand.

**The failures are not mostly vocabulary.** Measured over a two-day
session, they were: including working the reader cannot act on; asking for
a decision without supplying a recommendation; and using internal labels
(`H1`, `D2`) and self-coined names (*the probe*, *the walk*, *bake-off*)
as if they were shared. A term appearing in a **filename** is not shared
vocabulary. `reporting-in-chat` exists for exactly these.

The recurring failure is writing like an oracle — stating a conclusion in
terms the reader has to unpack, instead of showing the thing and saying
what it means. Three checks that catch it:

- **Never assert a property without its size, its subject, and its
  measure.** Not "the asymmetry is real and large" but "reads lose a
  repeat about five times as often as they gain one (2,438 against 501)".
  Words like *real*, *large*, *significant*, *agrees*, *tracks*,
  *dominates* are placeholders for a number — replace them.
- **Define a term the first time it does work in a sentence**, in the
  reader's language, before using it to carry an argument. If a sentence
  needs two undefined terms, it is two sentences.
- **When a claim can be checked against a figure, say where to look** —
  which panel, which bars, what the reader should compare. A claim about a
  plot that the plot appears to contradict is worse than no claim: state
  what is genuinely visible, then what it means.

Prefer natural frequencies to small percentages: "9 reads in 10,000"
lands where "0.09%" does not.

## Project layout pointers

- Architecture and specs: `ia/specs/` — start with
  `calling_pipeline_architecture.md` and `design_principles.md`.
- Implementation plans: `ia/feature_implementation_plans/`.
- Implementation reports: `ia/reports/implementations/` (create on first use).
- Skills used by the assistant: `ia/skills/`.
- Source tree: `src/` (library) and `src/main.rs` (CLI).
- Tests: per-module `#[cfg(test)] mod tests` blocks; integration
  tests in `tests/`.
