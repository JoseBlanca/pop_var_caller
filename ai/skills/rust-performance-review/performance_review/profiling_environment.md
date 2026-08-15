# Profiling environment: what runs where on this project's machines

This file is inventory, not method: which profiling tools are verified to work on this project's macOS host and its Linux dev container, what this hardware cannot measure, and how to diagnose a blocked sampling profiler on other machines. The orchestrator reads it at step 2 of the review procedure; every dispatched sub-agent reads it before writing a measurement plan, so that plans only name tools this environment can actually run.

Inventory entries are dated. If a run is far from the verification date, re-check with `command -v` — installs drift.

## Diagnosing a blocked sampling profiler (any machine)

**Linux** — `perf_event_open` is the kernel ABI all the major sampling tools use (`perf record`, `cargo flamegraph`, `samply`).

- `kernel.perf_event_paranoid >= 3` (Debian's default) blocks `perf_event_open` for unprivileged users. Fix: `sudo sysctl kernel.perf_event_paranoid=2` (persist via `/etc/sysctl.d/99-perf.conf`).
- Rootless container without `CAP_PERFMON`. Adding `--cap-add=PERFMON` on the container invocation only helps in *rootful* podman/docker — in rootless mode the kernel's `capable(CAP_PERFMON)` check runs against the host UID's caps, which the user namespace cannot grant. Lower the host paranoid instead, or run perf against the release binary from the host directly (`target/.../bench_binary --bench ... --profile-time 30`).
- Restricted CI / sandboxed runners. May require running benches on a privileged worker, or accepting that any code-level findings will be marked Speculative.

**macOS (incl. Apple Silicon M-series)** — There is no `perf_event_paranoid` knob; `perf` is Linux-specific. Profiling your own Rust process generally *just works* without elevated permissions. The native toolchain:

- `samply record cargo bench --bench <name> -- <filter> --profile-time 30` — same crate as Linux, cross-platform; outputs to a self-contained viewer (https://profiler.firefox.com). Easiest first try.
- `cargo flamegraph --bench <name> -- <filter> --profile-time 30` — on macOS this auto-switches to a `dtrace` backend.
- **Instruments.app** (ships with Xcode Command Line Tools, free) — the canonical GUI profiler; CPU Profiler template uses M-series hardware counters natively. CLI variant: `xcrun xctrace record --template 'Time Profiler' --launch -- <binary>`.
- `sample <pid>` — the simple, always-available CLI sampler that ships with the OS.

macOS blockers that *do* exist (less common, mention only if hit):

- System Integrity Protection (SIP) restricts some kernel-level DTrace probes; user-space sampling against your own process is unaffected. If `dtrace` complains about probe registration, disabling specific SIP features requires a Recovery-mode `csrutil` change — flag to user, don't recommend casually.
- Profiling a signed third-party binary (cross-process attach) needs the binary's entitlements to allow it. Not relevant for our own `cargo bench` output.

**Other OSes / WSL** — Treat as the closest match (WSL2 ≈ Linux; FreeBSD has `dtrace`). If you cannot match it, file the gap and ask the user.

## Tools verified available on this project's machine (macOS / Apple Silicon dev host)

Confirmed present as of 2026-06-01 (re-check with `command -v` if a run is far in the future — installs drift):

- **`samply`** (`~/.cargo/bin/samply`) — sampling profiler. `samply record --save-only -o out.json.gz -r 2000 -- <binary> …`. **Gotcha:** `--save-only` writes *unsymbolicated* module-relative addresses; symbolication happens in the Firefox-profiler UI (or via a symbol server samply runs). For a *headless* self-time ranking, prefer `sample` (below), which symbolicates inline.
- **`sample` (macOS built-in)** — attach to a running PID: `sample <pid> <secs> -file out.txt`. Produces a **symbolicated** call tree directly. This is the path that worked headlessly for the 2026-06-01 `var_calling` review. Note rayon pool threads show up as idle `Sleep::sleep` frames with the full sample count — filter them out.
- **`xctrace` (Instruments)** — `xcrun xctrace record --template 'Time Profiler' …`. GUI-oriented; `.trace` output is awkward to parse headlessly. Use `sample`/`samply` first.
- **host `cargo` + `rustc` (native arm64-apple-darwin)** — the dev container builds *Linux* binaries (`target-container/`), which cannot be profiled natively on the macOS host and whose in-VM `perf` targets the VM kernel. So for any host profiling you must **build natively on the host** (`cargo bench --no-run` / `cargo build --release --example …`), which lands in `target/` (gitignored) — this is the sanctioned exception to "cargo runs in the container". First native build is ~1 min.
- **DHAT** — no install needed; wired via the `dhat` crate behind the `dhat-heap` feature. `cargo run --release --example dhat_var_calling --features dhat-heap -- …` writes `dhat-heap.json` (open at <https://nnethercote.github.io/dh_view/dh_view.html>, or parse offline — the stacks are deep, so attribute to the first `src/<crate-module>/…rs:line` frame, skipping the `examples/dhat_*` alloc-hook frame and the alloc/core/BTreeMap internals). Sibling examples: `dhat_pileup`, `dhat_psp_reader`, `dhat_psp_writer`, `dhat_baq`. Run on the **host** when the fixture needs `$HOME/genomes/…` (the macOS container mounts only the project tree).
- **`cargo-show-asm`** (`cargo asm …`, installed 2026-06-01) — codegen inspection for the `hot_loops` category: confirm autovectorization / bounds-check elision on a named hot function (`cargo asm --lib --simplify "<crate>::path::to::fn"`). Functions that the profiler sees as distinct symbols are not fully inlined and `cargo asm` will find them; fully-inlined helpers must be inspected at their caller.
- **criterion benches** — `./scripts/dev.sh cargo bench --bench <name> -- <filter>` for reproducible wall numbers (run in the container for the committed baseline; run the native host build only when you also want to `samply`/`sample` it).
- **Not installed on the host** (it's in the container instead — see below): `cargo-flamegraph`, `perf`, `valgrind`, `coz`, `hyperfine`.

## Tools in the dev container (`./scripts/dev.sh` → Debian 12 / aarch64 — the production target)

Verified 2026-06-06. The container is **Linux aarch64** (the production build target), runs as **root**, `kernel.perf_event_paranoid=2`. **It is ephemeral** — every `./scripts/dev.sh …` is a fresh container; persist any built artifact or output under the bind-mounted project tree (`tmp/`), never `/tmp` (which is wiped). Already baked into the image (`Containerfile`): `perf` (`linux-perf`), `valgrind`, `hyperfine`, `cargo-flamegraph`, `samply`, `cargo-show-asm`. `coz` is in the `Containerfile` but **not yet in the live image** (its build layer lands on the next image rebuild — see note below); until then use the pre-built copy at `tmp/coz/install/bin/coz`.

> **PATH gotcha — use `bash -c`, not `bash -lc`.** A *login* shell sources `/etc/profile`, which resets `PATH` and **drops `/usr/local/cargo/bin`** — so `samply` / `cargo-flamegraph` / `cargo-asm` (all cargo-installed) silently vanish from a `bash -lc` command and appear "missing". `./scripts/dev.sh cargo …` and `./scripts/dev.sh bash -c '…'` both keep the right PATH. Reference genomes/PSPs are reachable via the same absolute host paths (HOME is propagated).
>
> **Image rebuilds:** `./scripts/rebuild-image.sh` (now runtime-aware: podman or Apple `container`). A *warm-cache* rebuild only re-runs changed layers (≈ minutes). A **cold** buildkit cache (fresh builder VM) forces a full from-scratch rebuild of this heavy image (samtools/bcftools/freebayes/GATK/nodejs + the Rust profilers + coz/LIEF) ≈ **40 min** — so don't trigger one casually; the coz/LIEF layer also bundles mbedTLS, which is the bulk of its build time.

- **`perf` CPU self-time sampling — WORKS, and on the production arch.** `perf record -e cpu-clock -F 499 -g -- <binary>` then `perf report`/`perf script`. Use the software `cpu-clock` event (not the default `cycles`).
- **`cargo flamegraph`** (already in the image) — wraps `perf`, so it needs the same `cpu-clock` override (its default `cycles` event returns `<not supported>` here). Working invocation (note `bash -c`, not `-lc`): `./scripts/dev.sh bash -c 'cargo flamegraph -c "record -e cpu-clock -F 997 --call-graph dwarf -g" --bin pop_var_caller -- var-calling --reference … --output /dev/null --threads 1 …'` → writes `flamegraph.svg`. The container is the right home for it (runs as root, no `sudo`; on the macOS host `cargo flamegraph` uses `dtrace`, which needs `sudo` and is awkward headless — prefer host `sample`/`samply` there instead). `samply` likewise wraps `perf` and takes the same `-e cpu-clock`.
- **`coz`** — **causal profiler**; the right tool for the producer→workers→writer pipeline ("speeding up which line actually moves end-to-end throughput?"). Add `#include <coz.h>` + a `COZ_PROGRESS` throughput point (e.g. per chunk written), build with `-g`, run `coz run --- <binary>`; reads `profile.coz`. **It samples on `PERF_TYPE_SOFTWARE`/`SW_TASK_CLOCK`, so it works here despite the missing PMU** (verified end-to-end). Until the next image rebuild bakes it in, run it from `tmp/coz/install/bin/coz` (set `COZ_LIBCOZ_PATH=tmp/coz/install/lib/aarch64-linux-gnu/libcoz.so` if invoked outside its prefix).
- **`valgrind`** — callgrind gives deterministic instruction-count ranking for A/B (slow, ~10–50×); cachegrind is unreliable on Apple Silicon — avoid.
- **`hyperfine`** — statistical CLI wall-clock benchmarking (reference-tool comparisons).

## What NONE of these environments can do (the PMU is not virtualized in the Apple-`container` VM)

`perf stat -e cache-misses,instructions,…` returns `<not supported>` in the container; consequently **`perf c2c` (false-sharing), `perf sched` / off-CPU blocking analysis, and hardware counters (IPC, cache/branch-miss) are unavailable** on this machine (host or container). This bites the concurrency category most (the "why doesn't wall scale past T≥2" question). Cover it by either:

1. **Instruments on the host** — the CPU Profiler / System Trace templates *do* read M-series hardware counters and thread-state (off-CPU), but are GUI-oriented and awkward to parse headlessly.
2. **`hotpath` software instrumentation** (https://hotpath.rs) — the right tool for *off-CPU* concurrency questions, and it needs **no PMU** (pure `Instant`/atomic instrumentation, not `perf_event_open`), so it works identically on host and in the Apple-`container` VM. Wrap crossbeam channels with `hotpath::channel!(…, wrap = true)` to get **current + max queue depth, max memory, and avg/p95 message delay** per channel; wrap locks with `hotpath::mutex!` / `hotpath::rw_lock!` to get **wait-before-acquire vs guard-hold time** (the "task holds the guard across a slow op; CPU profiler shows nothing, everyone else stalls" failure). This is the tool that finally *measures* the producer→workers→writer backpressure the cohort perf notes keep inferring indirectly (caller-starvation vs producer-back-pressure, "is the producer parked or contending?"). Wire it behind a `hotpath` cargo feature so the `hotpath::wrap::` types revert to the real std/crossbeam types and compile out of release builds. **Caveats:** `wrap = true` is crossbeam/std-only today and is **incompatible with `crossbeam::select!`** (needs a workaround if the stage uses it); wrapper types change channel signatures, so it is a feature-gated dev mode, not a drop-in. Complements `coz` (which answers *which line moving helps throughput*) rather than replacing it (`hotpath` answers *which stage is blocked and how deep its queue is*).
3. **Manual channel instrumentation** — time spent blocked in crossbeam `send` (workers can't keep up) vs `recv` (producer can't keep up); a few lines, always works, names the bottleneck stage directly. The zero-dependency fallback when you don't want to add the `hotpath` feature.
4. **A bare-metal Linux box** with full PMU + sched tracepoints if `perf c2c` / `perf sched` are truly needed.

**Branch-misprediction evidence without a branch-miss counter.** `perf stat -e branch-misses` is part of what the missing PMU takes away. The stand-ins are behavioural — run the same loop on the same values sorted vs shuffled, or sweep the predicate's selectivity (both described in `hot_loops.md`, branch rule) — plus `valgrind --tool=callgrind --branch-sim=yes`, which reports simulated conditional-branch and mispredict counts (`Bc`/`Bcm`). Callgrind simulates a simple predictor, not the M-series hardware, so treat its counts as comparative evidence between two variants of the same loop, not as absolute miss rates.

## Fixture locations

The reference FASTA for tomato fixtures lives at `~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`; real per-sample `.psp` cohorts are under `benchmarks/tomato1/results/ours/cohort/psp/` (and `/Users/jose/devel/pop_var_caller/tmp/aligned_psp/` for the 50-sample set).

## Toolchain pin

`rust-toolchain.toml` pins `channel = "1.95"` so criterion baselines stay comparable across review cycles — autovectorisation decisions shift between rustc versions without warning. Two consequences for measurement plans:

- A plan must not depend on a std API newer than 1.95 (e.g. `f64::algebraic_add`, stable from Rust 1.98) unless it names the pin bump as part of the fix's complexity cost.
- A pin bump is itself a measurable event: re-run the affected criterion baselines before and after, because codegen changes can move numbers with no code change.
