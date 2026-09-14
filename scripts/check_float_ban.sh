#!/bin/sh
# **Proves every entry of clippy.toml's float ban still refuses what it names.**
#
# `clippy.toml` refuses std's transcendental float methods so every such call goes through
# `crate::float` (doc/devel/implementation_plans/portable_float.md, step B6). Clippy silently
# ignores an entry whose path matches nothing — measured on Rust 1.98.0: a misspelled
# `f64::lnn` produced no warning at all — so a typo, or a method renamed upstream, would switch an
# entry off with every gate still green. Most entries never fire on the crate's own code, because
# nothing calls those methods; a silent entry and a working one look the same.
#
# This writes a throwaway example calling every banned method once, runs clippy on it, and fails
# unless clippy reports exactly the methods clippy.toml lists. The example is removed afterwards.
#
# Usage (run where cargo builds — in the dev container on a machine that has one):
#
#   ./scripts/dev.sh scripts/check_float_ban.sh
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
canary=examples/float_ban_canary.rs

# The banned paths, as clippy.toml lists them: `f64::ln`, `f32::powi`, ...
banned=$(grep -o 'path = "f[0-9]*::[a-z0-9_]*"' clippy.toml | sed 's/path = "//; s/"$//' | sort -u)
count=$(printf '%s\n' "$banned" | wc -l | tr -d ' ')

{
  echo "//! Throwaway canary written by scripts/check_float_ban.sh; never committed."
  echo "use std::hint::black_box;"
  echo "fn main() {"
  echo "    let x64: f64 = black_box(0.5);"
  echo "    let x32: f32 = black_box(0.5);"
  printf '%s\n' "$banned" | while IFS= read -r path; do
    ty=${path%%::*}
    method=${path#*::}
    case "$method" in
      powi) args="black_box(2)" ;;
      log | powf | atan2 | hypot) args="black_box(x$(echo "$ty" | tr -d f))" ;;
      *) args="" ;;
    esac
    echo "    let _ = black_box(x$(echo "$ty" | tr -d f)).$method($args);"
  done
  echo "}"
} > "$canary"
trap 'rm -f "$canary"' EXIT

output=$(cargo clippy --example float_ban_canary -- -D warnings 2>&1 || true)
# A misspelled or removed method does not compile, and then clippy reports nothing at all.
missing=$(printf '%s\n' "$output" | grep -o 'no method named `[a-z0-9_]*` found for type `f[0-9]*`' || true)
if [ -n "$missing" ]; then
  echo "float ban: clippy.toml names a method that does not exist (a typo or a rename):" >&2
  printf '%s\n' "$missing" >&2
  exit 1
fi
reported=$(printf '%s\n' "$output" | grep -o 'use of a disallowed method `f[0-9]*::[a-z0-9_]*`' \
  | sed 's/.*`\(.*\)`/\1/' | sort -u)
if [ "$reported" = "$banned" ]; then
  echo "float ban: all $count entries in clippy.toml refuse their method"
else
  echo "float ban: clippy.toml lists $count methods but clippy refused a different set" >&2
  echo "listed but not refused (a misspelled or renamed entry):" >&2
  printf '%s\n' "$banned" | grep -vxF "$reported" >&2 || true
  echo "refused but not listed:" >&2
  printf '%s\n' "$reported" | grep -vxF "$banned" >&2 || true
  printf '%s\n' "$output" | grep -E '^error(\[|:)' | grep -v 'disallowed method' | head -5 >&2
  exit 1
fi
