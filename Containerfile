# Development container for pop_var_caller.
#
# Serves two purposes:
#   1. Reproducible build/dev environment for the team.
#   2. Sandbox for running Claude Code with broad permissions, isolated
#      from the host filesystem outside the project directory.
#
# Rust version is pinned to match the host toolchain used during development.
FROM docker.io/library/rust:1.98-bookworm

# System dependencies:
#   - build-essential / pkg-config: required by a few cargo crates that
#     may compile native code.
#   - git / ca-certificates: needed by cargo for git-based deps and HTTPS.
#   - samtools / bcftools / tabix: inspect and compare BAM/VCF output
#     against reference tools.
#   - freebayes: reference variant caller; one of the three we compare
#     against in benchmarks/tomato1/.
#   - openjdk-17-jre-headless: Java runtime for GATK (installed below).
#   - unzip / wget: needed to fetch and unpack the GATK release.
#   - python3-psutil: lets the perf experiment scripts measure
#     subprocess RSS without needing uv inside the container.
#   - curl: bootstraps the NodeSource apt repo for Node.js.
#   - linux-perf: sampling profiler; backs cargo-flamegraph. In-container
#     sampling needs the host's kernel.perf_event_paranoid relaxed
#     (typically to 1). In rootless podman, --cap-add=SYS_ADMIN/PERFMON
#     does NOT help — the syscall still reaches the host kernel as the
#     unprivileged invoking user, so the sysctl is the only knob.
#   - hyperfine: statistical CLI benchmarking for reference comparisons.
#   - valgrind: callgrind/cachegrind for instruction-level profiling.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        git \
        ca-certificates \
        samtools \
        bcftools \
        tabix \
        freebayes \
        openjdk-17-jre-headless \
        unzip \
        wget \
        python3-psutil \
        curl \
        linux-perf \
        hyperfine \
        valgrind \
    && curl -fsSL https://deb.nodesource.com/setup_22.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# GATK 4: HaplotypeCaller + GenotypeGVCFs reference caller. Distributed
# as a Java fat-jar with a shell wrapper. Installs to /opt/gatk-<ver>/
# with a stable /opt/gatk -> /opt/gatk-<ver>/ symlink so the runner
# scripts (default `GATK_BIN=/opt/gatk/gatk`) stay valid across version
# bumps. Bump GATK_VERSION when you want to upgrade.
ARG GATK_VERSION=4.6.0.0
RUN wget -q "https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip" -O /tmp/gatk.zip \
    && unzip -q /tmp/gatk.zip -d /opt/ \
    && ln -s "/opt/gatk-${GATK_VERSION}" /opt/gatk \
    && rm /tmp/gatk.zip

# Rust components not included in the base image's minimal profile.
RUN rustup component add rustfmt clippy

# The three profiling tools, fetched as prebuilt binaries rather than compiled.
#
# **This was three `cargo install` layers and it was the bulk of a rebuild.**
# Each compiled its own dependency tree from source, and a rebuild that had to
# redo all three — which is what a base-image bump forces — ran for over an hour
# on the builder's default two cores. `cargo-binstall` pulls the upstream
# release binary where the project publishes one and falls back to compiling
# where it does not, so the slow path still exists but stops being the usual one.
#
# All three are pinned, as GATK and TRF-mod above are: an unpinned tool is a
# build that produces something different next month for no recorded reason.
# `--no-confirm` is required in a non-interactive build; `--locked` is kept for
# the fallback compile so a transitive bump cannot change what lands.
ARG BINSTALL_VERSION=1.10.17
ARG FLAMEGRAPH_VERSION=0.6.5
ARG SAMPLY_VERSION=0.13.1
ARG CARGO_SHOW_ASM_VERSION=0.2.47
RUN arch="$(uname -m)" \
    && case "$arch" in \
         aarch64) triple=aarch64-unknown-linux-musl ;; \
         x86_64)  triple=x86_64-unknown-linux-musl ;; \
         *) echo "no cargo-binstall build for $arch" >&2; exit 1 ;; \
       esac \
    && wget -q "https://github.com/cargo-bins/cargo-binstall/releases/download/v${BINSTALL_VERSION}/cargo-binstall-${triple}.tgz" \
         -O /tmp/binstall.tgz \
    && tar -xzf /tmp/binstall.tgz -C "$CARGO_HOME/bin" cargo-binstall \
    && rm /tmp/binstall.tgz \
    && cargo binstall --no-confirm --locked \
         "flamegraph@${FLAMEGRAPH_VERSION}" \
         "samply@${SAMPLY_VERSION}" \
         "cargo-show-asm@${CARGO_SHOW_ASM_VERSION}" \
    && flamegraph --version \
    && samply --version \
    && cargo asm --version

# coz: causal profiler — built for the cohort var-calling pipeline's
# producer->workers->writer shape, where "speed up which line actually
# improves end-to-end throughput?" is the question wall-clock and self-time
# profiles can't answer. Crucially it samples on PERF_TYPE_SOFTWARE /
# PERF_COUNT_SW_TASK_CLOCK, so it works inside the Apple-`container` Linux VM
# where the PMU is NOT virtualized (hardware counters, `perf c2c`, and
# `perf sched` off-CPU all return <not supported> there). Built from source
# pinned to a verified commit; CMake fetches LIEF + libelfin. Installs
# libcoz.so + the `coz` CLI to /usr/local. cmake is the only added build dep
# (build-essential / pkg-config are already present above).
#
# **Off by default since 2026-08-15, because it was most of the image's build
# time and nothing in the normal loop uses it.** coz answers one question — which
# line, if it got faster, would move end-to-end throughput — and it is the only
# layer here that compiles a large C++ project. Measured on a 10-core, 24 GB
# builder: everything up to this point finishes in about two and a half minutes,
# and this layer alone was still 30% through LIEF six minutes later. To get it:
#
#     ./scripts/rebuild-image.sh --build-arg WITH_COZ=1
#
# **The four-job cap is load-bearing whenever it *is* built.** LIEF's translation
# units want a couple of gigabytes of resident compiler each, so a bare `-j` —
# which takes the builder VM's core count — scales memory with cores and OOMs the
# moment the builder is given a real machine: on that same builder it died with
# `c++: fatal error: Killed signal terminated program cc1plus` on ELF/Builder.cpp,
# where the identical layer had completed on the old 2-core builder. Four holds at
# any builder size.
ARG WITH_COZ=0
RUN if [ "$WITH_COZ" != "1" ]; then \
        echo "coz: skipped (build with --build-arg WITH_COZ=1 to include it)"; \
    else \
        apt-get update \
        && apt-get install -y --no-install-recommends cmake \
        && git clone --recursive https://github.com/plasma-umass/coz /tmp/coz \
        && git -C /tmp/coz checkout 10630c542bc6d8a24595fd8283bea33bef892016 \
        && cmake -S /tmp/coz -B /tmp/coz/build -DCMAKE_BUILD_TYPE=Release \
        && cmake --build /tmp/coz/build -j4 \
        && cmake --install /tmp/coz/build \
        && ldconfig \
        && rm -rf /tmp/coz /var/lib/apt/lists/*; \
    fi

# TRF-mod (lh3): Tandem Repeats Finder with a BED-like output format — the
# SSR caller's Stage 0 (`ssr-catalog`) detection engine, shelled out per
# contig (see doc/devel/architecture/ssr_catalog.md §2). Single-file C build
# (`make -f compile.mak`, only build-essential needed), pinned to the commit
# the project vendors under TRF-mod/ so the catalog header's recorded
# `trf_mod_version` is reproducible. Installs to /usr/local/bin/trf-mod (on
# PATH for the catalog's layered binary discovery). Own layer so it does not
# invalidate the slow cargo-install / GATK layers above.
ARG TRF_MOD_COMMIT=3e891db310124f7e5f7a630a1c006650be9d1f3a
RUN git clone https://github.com/lh3/TRF-mod /tmp/trf-mod \
    && git -C /tmp/trf-mod checkout "${TRF_MOD_COMMIT}" \
    && make -C /tmp/trf-mod -f compile.mak \
    && install -m 0755 /tmp/trf-mod/trf-mod /usr/local/bin/trf-mod \
    && rm -rf /tmp/trf-mod

# Claude Code CLI, used when the container hosts an agent session.
RUN npm install -g @anthropic-ai/claude-code

# Pre-warm the cargo registry cache with this project's dependencies so that
# cold-start builds don't re-download every crate. Copy only the manifests
# (not the source) so this layer is cached until a dep changes. cargo fetch
# requires a parseable manifest plus src/main.rs, but it would also demand
# a file for every [[bench]] target — so we strip the [[bench]] sections
# from the manifest copy before fetching.
#
# **`[[example]]` sections are stripped for the identical reason, and were not
# until 2026-08-15.** The moment `Cargo.toml` gained an explicit `[[example]]`
# entry this layer stopped being buildable — `cargo fetch` refuses with
# "can't find `ng_generic_walk_probe` example at examples/…" because the
# `examples/` directory is not copied here either. A cached layer hid it, so it
# surfaced only on the next from-scratch build, which is the worst time to find
# it. Any target kind whose *files* this layer does not copy has to be stripped
# from the manifest it does. Only the populated registry under
# $CARGO_HOME survives into the final image; the stripped manifest, stub
# src/, and /build itself are byproducts of this layer.
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
RUN mkdir -p src \
    && echo 'fn main() {}' > src/main.rs \
    && awk 'BEGIN{skip=0} /^\[\[(bench|example)\]\]/{skip=1; next} /^\[/{skip=0} !skip' Cargo.toml > Cargo.toml.fetch \
    && mv Cargo.toml.fetch Cargo.toml \
    && cargo fetch \
    && rm -rf src

# GATK's `gatk` wrapper script uses `#!/usr/bin/env python`, but
# Debian ships only `python3` (no plain `python` symlink). Installed
# as its own layer so adding it doesn't invalidate the (slow) cargo
# install + GATK download layers above. If the Containerfile is ever
# rewritten, fold this back into the main apt-get block.
RUN apt-get update \
    && apt-get install -y --no-install-recommends python-is-python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY scripts/container-entrypoint.sh /usr/local/bin/container-entrypoint.sh
ENTRYPOINT ["/usr/local/bin/container-entrypoint.sh"]
CMD ["bash"]
