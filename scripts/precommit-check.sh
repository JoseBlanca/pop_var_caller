#!/usr/bin/env bash
# Pre-commit checks: fmt, clippy, tests, doc, bench-compile. Runs inside
# the dev container via scripts/dev.sh.
#
# Manual:       ./scripts/precommit-check.sh
# As git hook:  ln -s ../../scripts/precommit-check.sh .git/hooks/pre-commit
# Bypass once:  git commit --no-verify
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

"$SCRIPT_DIR/dev.sh" bash -c '
set -euo pipefail
step() { printf "\n\033[1;34m==> %s\033[0m\n" "$1"; }

step "1/5  cargo fmt --check"
cargo fmt --check

step "2/5  cargo clippy --all-targets --all-features -- -D warnings"
cargo clippy --all-targets --all-features -- -D warnings

step "3/5  cargo test"
cargo test

step "4/5  cargo doc  (a broken doc link is a broken promise)"
# **The same command and the same strictness as the doc job in CI**, so that a
# doc link this catches is one CI would have caught and not a stricter local
# rule. Under RUSTDOCFLAGS=-D warnings the warnings from rustdoc count too: a
# link whose explicit target is redundant fails here, which is what keeps the
# habit of writing one from spreading.
#
# ⚠ This whole block is one single-quoted argument to bash -c, so an apostrophe
# anywhere in it ends the string. Write around them.
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib --all-features

step "5/5  cargo bench --no-run  (compile-only, catches bench bitrot)"
cargo bench --no-run

printf "\n\033[1;32mAll checks passed.\033[0m\n"
'
