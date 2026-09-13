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

step "1/6  dependency guard: production must not depend on ng"
# ng is a from-scratch experiment caller; production (everything outside
# src/ng/) must never `use crate::ng`. The exp binary
# (src/pop_var_caller_exp/) is a driver for ng, so it is allowed to depend
# on ng and is excluded too — the dependency points the safe way (see
# doc/devel/ng/spec/typed_regions_cli.md T7). The excluded paths are
# load-bearing: widen this list and the guard becomes a comment.
guard_hits=$(grep -rn "use crate::ng" src --include="*.rs" \
    | grep -Ev "^src/(ng|pop_var_caller_exp)/" || true)
if [ -n "$guard_hits" ]; then
    printf "\033[1;31mguard failed: production code depends on ng:\033[0m\n%s\n" "$guard_hits"
    exit 1
fi

step "1b/6 dependency guard: shipped ng must not depend on production"
# **The mirror of the step above, added by the promotion**
# (doc/devel/implementation_plans/promote_ng_to_production.md, milestone B). Nothing
# under src/ng/ that a run executes may use a module production owns, because those
# modules are being deleted. What ng needed from them has been copied in.
#
# **Two rules together decide what counts as shipped**, because neither alone does.
#
# Column zero: a module-level "use" starts there, while a use inside a test module or
# a test function is indented. That excludes the oracle imports embedded in ordinary
# modules, of which there are about twenty.
#
# The file list: a whole module gated #[cfg(test)] at its parent has its imports at
# column zero like any other, so the dedicated parity and copy-fidelity files need
# naming. Each one exists to run production beside ng and assert the two agree, which
# is the whole reason the promotion severs before it deletes; they keep reaching in
# until milestone C freezes their answers into fixtures, and the list empties then.
#
# Matching on indentation is a heuristic and its limit is worth stating: a shipped
# import written indented inside a function would slip past. rustfmt does not produce
# one at module level, and the compiler is the real check the moment production is
# gone. This exists so the end state of milestone B is visible before D starts, and
# so that a new reach-in fails here rather than at the deletion.
production_modules="pileup|psp|pileup_record|genetics|pop_var_caller|var_calling|vcf|ssr|paralog|sample_summary|baq|norm_seqs"
oracles="parity\.rs|leftmost_property\.rs|test_fixtures\.rs|mock_reference\.rs|ssr_production_differential\.rs|prepared_read\.rs"
reach_in=$(grep -rnE "^use crate::($production_modules)(::|;| )" src/ng --include="*.rs" \
    | grep -Ev "$oracles" || true)
if [ -n "$reach_in" ]; then
    printf "\033[1;31mguard failed: shipped ng code depends on production:\033[0m\n%s\n" "$reach_in"
    printf "copy what is needed into ng, as milestone B did, rather than reaching in.\n"
    exit 1
fi

step "2/6  cargo fmt --check"
cargo fmt --check

step "3/6  cargo clippy --all-targets --all-features -- -D warnings"
cargo clippy --all-targets --all-features -- -D warnings

step "4/6  cargo test"
cargo test

step "5/6  cargo doc  (a broken doc link is a broken promise)"
# **The same command and the same strictness as the doc job in CI**, so that a
# doc link this catches is one CI would have caught and not a stricter local
# rule. Under RUSTDOCFLAGS=-D warnings the warnings from rustdoc count too: a
# link whose explicit target is redundant fails here, which is what keeps the
# habit of writing one from spreading.
#
# ⚠ This whole block is one single-quoted argument to bash -c, so an apostrophe
# anywhere in it ends the string. Write around them.
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib --all-features

step "6/6  cargo bench --no-run  (compile-only, catches bench bitrot)"
cargo bench --no-run

printf "\n\033[1;32mAll checks passed.\033[0m\n"
'
